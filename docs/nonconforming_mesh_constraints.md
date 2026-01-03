# 非共形网格的约束与变换矩阵机制

## 1. 概述

本文档详细说明在非共形网格（如mortar方法、挂节点AMR）中如何在刚度矩阵和线性方程组上添加约束和变换矩阵。这与共形多patch NURBS的DOF共享机制形成对比。

## 2. 非共形网格的类型

### 2.1 Mortar方法（不匹配网格）
- 两个网格在界面上不一致（节点位置不重合）
- 需要通过投影或插值建立DOF关系

### 2.2 挂节点AMR
- 局部细化导致相邻单元大小不同
- 细化边上的节点是"挂节点"，需要约束

### 2.3 滑动界面
- 两个网格区域可相对滑动
- 界面关系随时间变化

## 3. 约束方程的数学形式

### 3.1 基本约束
设细边上的挂节点DOF为 $u_s$（slave），粗边上的DOF为 $u_m$（master），约束关系为：

$$u_s = C \cdot u_m$$

其中 $C$ 是约束矩阵（或变换矩阵），通常通过插值或L2投影获得。

### 3.2 示例：2D挂节点

```
粗边:     o---------o
          1         2

细边:     o----o----o
          1   (s)   2
```

挂节点 $s$ 的约束：
$$u_s = 0.5 \cdot u_1 + 0.5 \cdot u_2$$

约束矩阵：$C = [0.5, 0.5]$

## 4. 刚度矩阵的修改

### 4.1 原始系统

$$K \cdot u = f$$

其中 $K$ 是刚度矩阵，$u$ 是所有DOF（包括slave和master），$f$ 是载荷向量。

### 4.2 引入约束后的分块形式

将DOF分为独立DOF $u_m$ 和从属DOF $u_s$：

$$\begin{bmatrix} K_{mm} & K_{ms} \\ K_{sm} & K_{ss} \end{bmatrix} \begin{bmatrix} u_m \\ u_s \end{bmatrix} = \begin{bmatrix} f_m \\ f_s \end{bmatrix}$$

约束：$u_s = C \cdot u_m$

### 4.3 方法一：约束消除（静态凝聚）

用约束替换 $u_s$，得到缩减系统：

$$\tilde{K} \cdot u_m = \tilde{f}$$

其中：
- $\tilde{K} = K_{mm} + K_{ms} C + C^T K_{sm} + C^T K_{ss} C$
- $\tilde{f} = f_m + C^T f_s$

**实现代码框架：**
```cpp
// 构建缩减刚度矩阵
SparseMatrix K_reduced;
K_reduced = K_mm;
K_reduced.Add(1.0, K_ms * C);
K_reduced.Add(1.0, Ct * K_sm);
K_reduced.Add(1.0, Ct * K_ss * C);

// 缩减载荷向量
Vector f_reduced = f_m;
Ct.AddMult(f_s, f_reduced);

// 求解缩减系统
Solve(K_reduced, u_m, f_reduced);

// 恢复从属DOF
C.Mult(u_m, u_s);
```

### 4.4 方法二：Lagrange乘子法

增广系统，引入Lagrange乘子 $\lambda$ 强制约束：

$$\begin{bmatrix} K & G^T \\ G & 0 \end{bmatrix} \begin{bmatrix} u \\ \lambda \end{bmatrix} = \begin{bmatrix} f \\ 0 \end{bmatrix}$$

其中 $G = [-C \quad I]$ 使得 $G \cdot u = u_s - C u_m = 0$

**优点**：系统对称，可使用标准求解器
**缺点**：系统规模增大，需要合适的预条件

### 4.5 方法三：惩罚法

添加惩罚项近似强制约束：

$$K_\text{aug} = K + \beta G^T G$$

其中 $\beta$ 是大的惩罚参数。

**优点**：不增加系统规模
**缺点**：病态性，$\beta$ 选择困难

## 5. MFEM中的实现

### 5.1 非共形AMR（NCMesh）

MFEM使用约束矩阵 `P` 来处理挂节点：

```cpp
// fem/fespace.cpp
if (mesh->ncmesh)
{
   // P 是插值矩阵，从共形DOF到所有DOF
   P = GetConformingProlongation();
   R = GetConformingRestriction();
   
   // 实际求解在共形空间进行
   // u_all = P * u_conforming
}
```

### 5.2 约束矩阵的组装

```cpp
// mesh/ncmesh.cpp
void NCMesh::GetFineToCoarseConstraints(
    int entity_type,
    Table &constraints) const
{
    // 为每个挂节点计算其与粗节点的约束系数
    for (auto& hanging_node : hanging_nodes)
    {
        // 找到父边/面的节点
        Array<int> master_dofs;
        GetMasterDofs(hanging_node, master_dofs);
        
        // 计算插值系数
        Vector coeffs;
        ComputeInterpolationCoeffs(hanging_node, master_dofs, coeffs);
        
        constraints.AddRow(hanging_node, master_dofs, coeffs);
    }
}
```

### 5.3 RAP变换

对于稀疏矩阵，使用RAP（Restriction-Apply-Prolongation）变换：

```cpp
// A 是原始矩阵（包含挂节点）
// P 是prolongation矩阵
// A_c 是缩减后的共形矩阵

SparseMatrix *A_c = RAP(P, A, P);  // P^T * A * P

// 求解
Solve(A_c, x_c, b_c);

// 恢复完整解
P->Mult(x_c, x);
```

## 6. H(curl)空间的特殊处理

### 6.1 切向连续性约束

H(curl)空间要求切向分量连续，约束更复杂：

```cpp
// 对于H(curl)，约束需要考虑边的方向
// 如果边方向相反，约束系数需要取负

if (edge_orientation < 0)
{
    constraint_coeff *= -1.0;
}
```

### 6.2 向量插值

挂节点处的向量场插值：

$$\mathbf{E}_s = \sum_i c_i \mathbf{E}_{m_i}$$

需要确保切向分量正确插值。

## 7. 共形vs非共形对比

| 特性 | 共形多Patch NURBS | 非共形网格 |
|------|-------------------|------------|
| DOF处理 | 直接共享 | 约束关系 |
| 矩阵组装 | 标准累加 | 需要变换矩阵 |
| 系统规模 | 不变 | 可能增大（Lagrange） |
| 数值精度 | 精确匹配 | 依赖插值精度 |
| 实现复杂度 | 简单 | 复杂 |

## 8. 代码示例：非共形约束的完整流程

```cpp
// 1. 检测非共形边界
NCMesh *ncmesh = mesh->ncmesh;

// 2. 构建约束矩阵
SparseMatrix *P = fespace->GetConformingProlongation();
SparseMatrix *R = fespace->GetConformingRestriction();

// 3. 组装原始刚度矩阵（包含所有DOF）
BilinearForm a(fespace);
a.AddDomainIntegrator(new DiffusionIntegrator);
a.Assemble();
SparseMatrix A_full = a.SpMat();

// 4. 应用RAP变换得到共形系统
SparseMatrix *A_conf = RAP(*P, A_full, *P);

// 5. 变换右端项
Vector b_full, b_conf;
// ... 组装 b_full ...
R->Mult(b_full, b_conf);

// 6. 求解共形系统
Vector x_conf;
PCG(*A_conf, b_conf, x_conf);

// 7. 恢复完整解
Vector x_full;
P->Mult(x_conf, x_full);
```

## 9. 总结

### 9.1 非共形网格的约束机制

1. **约束方程**：$u_s = C \cdot u_m$ 建立slave和master DOF的关系
2. **刚度矩阵变换**：通过 $\tilde{K} = P^T K P$ 或Lagrange乘子法
3. **线性方程组**：在缩减空间求解，再恢复完整解

### 9.2 与共形多Patch的区别

共形多patch NURBS不需要这些机制，因为：
- 共享边上的DOF天然相同（相同全局编号）
- 刚度矩阵组装时贡献直接累加
- 无需约束矩阵或变换

### 9.3 选择建议

- **共形网格**：优先使用，实现简单，精度高
- **非共形网格**：当几何或网格约束需要时使用，需要额外的约束处理
