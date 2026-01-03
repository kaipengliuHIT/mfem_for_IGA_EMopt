# 多Patch H(curl) NURBS 修复文档

## 1. Bug 原因分析

### 1.1 问题现象
在多patch NURBS网格上使用H(curl)有限元空间时，线性求解器无法收敛，报告"矩阵不正定"错误。

### 1.2 根本原因

#### 问题定位
通过诊断发现：
- **DOF表大小与有限元DOF数量不匹配**：例如FE期望24个DOF，但DOF表有28个条目
- **矩阵对角线有零值**：某些DOF的对角线元素为零或接近零

#### 根本原因
`GetCurlExtension()`函数中的order计算逻辑错误。

**原始代码：**
```cpp
Array<int> newOrders = GetOrders();
for (int c = 0; c < newOrders.Size(); c++) { newOrders[c]++; }
newOrders[component] -= 1;
```

**问题分析：**
- `GetOrders()` 返回所有唯一KnotVector的orders
- 在多patch情况下，唯一KnotVector数量可能大于维度
  - 2D两patch: 可能有3个唯一KnotVector（共享边合并）
  - 3D两patch: 可能有4-6个唯一KnotVector
- 原代码简单地对所有KnotVector应用`+1`再对`component`位置`-1`
- 这导致了错误的order分配，因为KnotVector索引并不直接对应方向

**示例（2D两patch，原代码）：**
```
原始orders: [2, 2, 2]  (3个唯一KV)
错误计算 component=0: [2, 3, 3]  // 只有index 0保持不变
正确应该: [2, 3, 2]  // 所有x方向KV保持2，y方向KV加1
```

### 1.3 影响
- Component 0 (x方向): orders=(2,3,3,2) 应该是 (2,3,2,3)
- DOF表每个元素的DOF数量错误（16而不是12）
- 合并后总DOF数与FE期望不匹配
- 某些DOF的形函数未被正确计算，导致矩阵对角线为零

## 2. 修复内容

### 2.1 核心修复：`GetCurlExtension()` 函数

**文件：** `/mesh/nurbs.cpp`

**修复思路：**
使用`edge_to_ukv`映射来确定每个唯一KnotVector对应的方向，然后正确应用H(curl)的order变换。

**修复代码：**
```cpp
NURBSExtension* NURBSExtension::GetCurlExtension(int component)
{
   const int dim = Dimension();
   const int np = GetNP();
   Array<int> newOrders = GetOrders();
   
   // 构建唯一KnotVector索引到方向的映射
   Array<int> ukv_to_dir(newOrders.Size());
   ukv_to_dir = -1;
   
   // 遍历所有patch的边，确定每个KnotVector的方向
   for (int p = 0; p < np; p++)
   {
      Array<int> edges, orient;
      patchTopo->GetElementEdges(p, edges, orient);
      
      for (int e = 0; e < edges.Size(); e++)
      {
         int ukv_idx = edge_to_ukv[edges[e]];
         int dir;
         if (dim == 2)
         {
            // 2D: edges 0,2 是 x方向, edges 1,3 是 y方向
            dir = (e == 0 || e == 2) ? 0 : 1;
         }
         else // dim == 3
         {
            // 3D: edges 0,2,4,6 是 x, 1,3,5,7 是 y, 8-11 是 z
            if (e < 4) { dir = e % 2; }
            else if (e < 8) { dir = (e - 4) % 2; }
            else { dir = 2; }
         }
         
         if (ukv_to_dir[ukv_idx] == -1)
         {
            ukv_to_dir[ukv_idx] = dir;
         }
      }
   }
   
   // 基于方向应用H(curl) order变换
   for (int ukv = 0; ukv < newOrders.Size(); ukv++)
   {
      int dir = ukv_to_dir[ukv];
      if (dir != component)
      {
         newOrders[ukv]++;  // 非component方向 +1
      }
   }

   return new NURBSExtension(this, newOrders, Mode::H_CURL);
}
```

### 2.2 其他相关修改

1. **移除单patch限制**：删除了`GetCurlExtension`中的单patch检查
2. **DOF连接框架**：实现了`AutoConnectPatchBoundariesHCurl`和`ConnectPatchEdgeHCurl2D`（虽然最终发现mesh拓扑已处理DOF共享）
3. **测试用例**：`nurbs_ex25p_multipatch.cpp`用于测试多patch H(curl)

## 3. 兼容性说明

### 3.1 单进程运行
✅ **完全兼容**

```bash
mpirun -np 1 ./nurbs_ex25p_multipatch -m square-nurbs.mesh
# 结果: L2误差 = 0.00337522
```

### 3.2 单patch网格
✅ **完全兼容**

修复后的代码对单patch网格完全透明：
- 单patch时，`ukv_to_dir`映射正确退化为简单的方向索引
- 所有单patch测试通过

### 3.3 多patch 2D
✅ **已修复并测试**

```bash
mpirun -np 2 ./nurbs_ex25p_multipatch -m two-squares-nurbs.mesh
# 结果: L2误差 = 0.00477329
```

### 3.4 多patch 3D
✅ **已修复并测试**

```bash
mpirun -np 1 ./nurbs_ex25p_multipatch -m two-cubes-nurbs.mesh
# 结果: L2误差 = 0.00580077
```

### 3.5 三个及更多patch
⚠️ **理论支持，需要测试**

修复逻辑是通用的：
- 遍历所有patch的边来构建方向映射
- 不依赖于patch数量
- 只要mesh拓扑正确（共享边/面正确识别），应该能正常工作

**建议测试：**
- 使用`plus-nurbs.mesh`（多于2个patch的2D网格）
- 创建三个cube的3D网格进行测试

## 4. 技术细节

### 4.1 KnotVector组织

MFEM中有两种KnotVector存储：

1. **`knotVectors`**: 唯一KnotVector数组（共享的KV只存储一次）
2. **`knotVectorsCompr`**: 综合KnotVector数组，按`[patch0_dir0, patch0_dir1, ..., patch1_dir0, ...]`组织

`GetOrders()`返回`knotVectors`的orders，而`edge_to_ukv`提供从边索引到唯一KV索引的映射。

### 4.2 边到方向的映射规则

**2D (四边形patch):**
- 边 0, 2: x方向
- 边 1, 3: y方向

**3D (六面体patch):**
- 边 0,2,4,6: x方向
- 边 1,3,5,7: y方向  
- 边 8,9,10,11: z方向

### 4.3 H(curl) Order规则

对于H(curl)空间的component方向`c`：
- 方向`c`的order保持不变
- 其他方向的order +1

例如2D order 2的H(curl)：
- Component 0 (x): orders = (2, 3)
- Component 1 (y): orders = (3, 2)

## 5. 测试命令

```bash
# 构建
cd /path/to/mfem_hcurl_multipatch/build
make -j4

# 构建测试
cd /path/to/mfem_hcurl_multipatch/miniapps/nurbs
make nurbs_ex25p_multipatch

# 测试单patch
mpirun -np 1 ./nurbs_ex25p_multipatch -m ../../data/square-nurbs.mesh

# 测试2D多patch
mpirun -np 2 ./nurbs_ex25p_multipatch -m meshes/two-squares-nurbs.mesh

# 测试3D多patch
mpirun -np 1 ./nurbs_ex25p_multipatch -m meshes/two-cubes-nurbs.mesh
```

## 6. 总结

| 场景 | 状态 | L2误差 |
|------|------|--------|
| 单patch 2D | ✅ | 0.00337522 |
| 多patch 2D (2 patches) | ✅ | 0.00477329 |
| 多patch 3D (2 patches) | ✅ | 0.00580077 |
| 并行运行 (np=2) | ✅ | 0.00477329 |
| 三个以上patch | ⚠️ 需测试 | - |

修复的关键是正确建立唯一KnotVector索引到空间方向的映射，而不是简单地使用索引模维度。
