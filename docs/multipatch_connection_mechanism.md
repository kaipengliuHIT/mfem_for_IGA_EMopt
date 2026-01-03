# 多Patch NURBS 连接机制详解

## 1. 概述

在MFEM的多patch NURBS实现中，patch之间的连接**不是**通过在刚度矩阵上添加变换矩阵（如非共形网格的约束方程），而是通过**DOF共享机制**实现的。

## 2. DOF共享机制

### 2.1 基本原理

多patch NURBS网格中，相邻patch共享边界上的控制点。MFEM通过以下方式处理：

1. **mesh拓扑层面**：共享边/面上的顶点被识别为同一顶点
2. **NURBSExtension层面**：共享边上的DOF被映射到相同的全局DOF编号
3. **刚度矩阵层面**：共享DOF只出现一次，贡献从两个patch累加

### 2.2 与非共形网格的区别

| 特性 | 多Patch NURBS | 非共形网格 |
|------|---------------|------------|
| 连接方式 | DOF共享（共形） | 约束方程/变换矩阵 |
| 刚度矩阵 | 直接累加贡献 | 需要额外约束项 |
| DOF数量 | 共享DOF只计一次 | 每侧独立DOF |
| 连续性 | 自然满足 | 通过约束强制 |

## 3. 刚度矩阵组装过程

### 3.1 元素循环

```
对于每个元素 e:
    1. 获取元素DOF表: fespace->GetElementDofs(e, dofs)
    2. 计算局部刚度矩阵: K_e
    3. 累加到全局矩阵: A(dofs[i], dofs[j]) += K_e(i,j)
```

### 3.2 共享DOF的处理

当两个patch共享边界时：
- **Patch 0 的边界元素**：其DOF表包含共享边上的DOF编号（如 DOF 5, 6, 7）
- **Patch 1 的边界元素**：其DOF表包含**相同的**DOF编号（DOF 5, 6, 7）

因此，在矩阵组装时：
```
A(5,5) = K_patch0(i,i) + K_patch1(j,j)  // 两个patch的贡献自动累加
```

### 3.3 示例：两个正方形patch

```
Patch 0 (左)          Patch 1 (右)
+---+---+            +---+---+
| 0 | 1 |            | 2 | 3 |
+---+---+            +---+---+
| 4 | 5 |            | 6 | 7 |
+---+---+            +---+---+
    共享边
```

共享边上的DOF在两个patch中使用相同的全局编号，例如：
- Patch 0 元素 1 的右边界DOF: [10, 11, 12]
- Patch 1 元素 2 的左边界DOF: [10, 11, 12]  ← 相同！

## 4. MFEM实现细节

### 4.1 关键数据结构

```cpp
class NURBSExtension {
    Array<int> edge_to_ukv;      // 边到唯一KnotVector的映射
    Array<int> d_to_d;           // DOF到DOF的映射（用于周期性BC）
    Array<KnotVector*> knotVectors;      // 唯一KnotVector集合
    Array<KnotVector*> knotVectorsCompr; // 综合KnotVector（按patch组织）
};
```

### 4.2 DOF编号过程

1. **GenerateElementDofTable()**：为每个元素生成局部DOF表
2. **NURBSPatchMap::operator()**：将patch局部(i,j)索引映射到全局DOF
3. **共享边处理**：通过`v_spaceOffsets`和`e_spaceOffsets`确保共享实体有相同DOF

### 4.3 NURBSPatchMap 映射逻辑

```cpp
// 2D情况下的DOF映射
int NURBSPatchMap::operator()(const int i, const int j) const
{
    switch (3*F(j, J) + F(i, I))
    {
        case 0: return verts[0];      // 左下角顶点
        case 1: return EC(0, i, I);   // 底边
        case 2: return verts[1];      // 右下角顶点
        case 3: return EC(3, j, J);   // 左边
        case 4: return pOffset + ...  // 内部DOF
        case 5: return EC(1, j, J);   // 右边
        case 6: return verts[3];      // 左上角顶点
        case 7: return EC(2, i, I);   // 顶边
        case 8: return verts[2];      // 右上角顶点
    }
}
```

关键是：`EC()`函数返回的边DOF编号，对于共享边，两个patch会返回相同的值。

## 5. H(curl) 空间的特殊处理

### 5.1 向量空间的DOF表合并

H(curl)空间由多个标量组件组成：
```cpp
// fespace.cpp 中的合并逻辑
elem_dof = new Table(*VNURBSext[0]->GetElementDofTable(),
                     *VNURBSext[1]->GetElementDofTable(), offset1);
```

每个组件的DOF表被拼接，第二个组件的DOF编号加上offset。

### 5.2 共享边上的H(curl) DOF

对于H(curl)，只有**切向分量**需要在共享边上连续：
- 水平边（y=常数）：x分量（Component 0）需要连续
- 垂直边（x=常数）：y分量（Component 1）需要连续

这通过每个组件独立的NURBSExtension自动处理。

## 6. 验证共享DOF

### 6.1 检查DOF表

```cpp
// 打印两个相邻元素的DOF表
Array<int> dofs0, dofs1;
fespace->GetElementDofs(elem_patch0, dofs0);
fespace->GetElementDofs(elem_patch1, dofs1);

// 找出共享的DOF（应该有重叠）
for (int i = 0; i < dofs0.Size(); i++) {
    for (int j = 0; j < dofs1.Size(); j++) {
        if (dofs0[i] == dofs1[j]) {
            cout << "Shared DOF: " << dofs0[i] << endl;
        }
    }
}
```

### 6.2 检查刚度矩阵

共享DOF的对角线元素应该包含来自两个patch的贡献：
```cpp
// 对于共享DOF d，其对角线值应该是两个patch贡献之和
A(d,d) = A_patch0(d,d) + A_patch1(d,d)
```

## 7. 测试结果

| 网格 | Patch数 | DOF数 | 收敛 | L2误差 |
|------|---------|-------|------|--------|
| two-squares-nurbs.mesh | 2 | 71 | ✅ | 0.00477329 |
| plus-nurbs.mesh | 5 | 164 | ✅ | 0.00754723 |
| two-cubes-nurbs.mesh | 2 (3D) | 535 | ✅ | 0.00580077 |

## 8. 总结

多patch NURBS的连接机制是**共形的DOF共享**，不需要额外的变换矩阵或约束方程。这与非共形网格（如mortar方法）有本质区别：

1. **优点**：自然满足连续性，无需额外约束
2. **要求**：patch边界必须几何上重合，且阶数一致
3. **实现**：通过NURBSPatchMap确保共享边/面的DOF编号一致

刚度矩阵的组装过程与单patch完全相同，共享DOF的贡献自动累加，无需特殊处理。
