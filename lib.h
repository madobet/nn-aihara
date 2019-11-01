/*!
 * @file lib.h
 * @brief 各种库函数
 * @details 目前包括 IO类、Math类和Graph类
 * @mainpage Library
 * @author Tooko Madobe
 * @email madobet@outlook.com
 * @version 0.1
 * @date 2019-03-29
 * @license GPLv3
 */
#ifndef NTOOLKIT_LIB_H
#define NTOOLKIT_LIB_H

#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <random>

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <igraph/igraph.h>

namespace ntoolkit{

  /*! @brief 输入输出方法类
   *  处理程序的各种各样输入输出
   */
  class IO{
  public:
      /*!
       * @brief 将 str 中的 old_str 字符串全部用 new_str 替换掉
       * @details 搜寻到字符串末尾跳出循环
       * @param str 
       * @param old_str
       * @param new_str
       * @return 修改后原字符串的引用（注意此时原字符串已经被修改）
       *    @retval std::string&
       */
      static std::string& ReplaceAllStr(std::string& str, const std::string& old_str, const std::string& new_str);

      /*!
       * @brief 从文件自动分析矩阵维数并读入矩阵
       * @tparam T
       * @param filepath
       * @return 读入的矩阵
       *    @retval std::vector< std::vector<T> >
       */
      template <typename T>
      static std::vector< std::vector<T> > ReadMatFromFile(const std::string &filepath);

      /*!
       * @brief 从文件自动分析矩阵维数并以一定阈值读入矩阵
       * @tparam T
       * @param filepath
       * @param threshold
       * @return 读入的矩阵
       *    @retval std::vector< std::vector<T> >
       */
      template <typename T>
      static std::vector< std::vector<T> > ReadMatFromFile(const std::string &filepath, T threshold);

      /*!
       * @brief 从矩阵文件取得矩阵维数
       * @param filepath 文件路径
       * @return 矩阵的维数
       *    @retval int
       */
      static int GetMatDimFromFile(const std::string &filepath);

      /*!
       * @brief 将向量写入指定文件
       * @tparam T 向量的数据类型 int 或 double
       * @param array 向量
       * @param filename 文件名
       * @return 写入成功或失败
       *    @retval bool
       */
      template <typename T>
      static bool SaveVec2File(std::vector<T> &array, std::string &filename);

      /*!
       * @brief 将矩阵写入指定文件
       * @tparam T 矩阵的数据类型 int 或 double
       * @param m 矩阵
       * @param fpath 文件路径
       * @return 写入成功或失败
       *    @retval bool
       */
      template <typename T> static bool WriteMat2File(std::vector<std::vector<T> > &m, std::string &fpath);

      /*!
       * @brief 从矩阵中取出一列
       * @tparam T 矩阵的数据类型 int 或 double
       * @param m 矩阵
       * @param cid 列数
       * @return (列)向量
       *    @retval std::vector<T>
       */
      template <typename T> std::vector<T> ExtOneColFromMat(std::vector<std::vector<T> > &m, int cid);
  };

  /*!
   * @brief 数学库
   * 数学常量、函数和算法
   */
  class Math{
  public:
      static const int max_buffer_size = 512;
      static const int infinite_int = 0x7FFF;
      static const double pi;
      static const double eps;
      template <typename T>
      static inline T Abs(T x) { return (x<0)?(-x):x;}
      template <typename T>
      static inline T Sgn(T x) { return (x>0)?1:( (x<0)?(-1):0 ); }
      static inline double Deg2Rad(double deg) { return (2*pi*deg)/360; }
      static inline double Sigmoid(double x) { return 1/(1 + exp(-x/eps)); }

      /*!
       * @brief 生成指定范围(pair)内整型随机数
       * @details 不采用模板由于 uniform_int_distribution uniform_real_distribution 名字不一致
       * @param lim_lo2up
       * @return 整型随机数
       *    @retval int
       */
      static int RandNum(std::pair<int, int> &lim_lo2up);

      /*!
       * @brief 生成指定范围(pair)内浮点随机数
       * @details 不采用模板由于 uniform_int_distribution uniform_real_distribution 名字不一致
       * @param lim_lo2up 范围
       * @return 浮点随机数
       *    @retval double
       */
      static double RandNum(std::pair<double, double> &lim_lo2up);

      /*!
       * @brief 生成随机向量
       * @tparam T
       * @param size 向量维数
       * @param lim_lo2up 范围
       * @return size 维的随机值向量
       *    @retval std::vector<T>
       */
      template <typename T>
      static std::vector<T> RandVec(T size, std::pair<T, T> &lim_lo2up);
      
      /*!
       * @brief 生成m*n的随机矩阵
       * @tparam T
       * @param m_mul_n m和n
       * @param lim_lo2up 范围
       * @return m*n具有随机值的矩阵
       *    @retval std::vector< std::vector<T> >
       */
      template <typename T>
      static std::vector< std::vector<T> > RandMat(std::pair<T, T> &m_mul_n, std::pair<T, T> &lim_lo2up);

      /*!
       * @brief 简易计算器（CalFunc的后端）
       * @details 计算算数表达式的值
       *
       * Usage
       * - support: + - * / sin cos tan sqrt pow log ln exp
       * - must using operator!
       * - operator:+ - * / s   c   t   v    ^   L   l  e
       * - 小数运算，运算数据类型为double（双精度）型
       * - 小括号（）
       *
       * example:
       * - a^b 表示求a的b次方
       * - va或av 表示求a平方根
       * - v(va)与avv等价 表示对a求两次平方根，可以此类推，求a的偶次方根。
       * - s30 表示求30°的正弦值
       * - L100 表示求100以10为底的对数, L要大写
       *
       * @warning 不支持负数：如果要求表达式 -5*9 ,则可以写成 (0-5)*9 。总之可以用括号解决
       * @warning 不进行表达式的正确性验证
       *
       * 算术优先级:
       * log~ln > sqrt~^ > sin~cos~tan > *~/ > +~-
       * 
       * @param expression 表达式
       * @return 计算值
       *    @retval double
       */
      double CalExp(char *expression);

      /*!
       * @brief 函数计算器
       * @param func 函数表达式
       * @param arg 自变量字母
       * @param arg_value 自变量值
       * @return 函数的计算值
       *    @retval double
       * @warning arg自变量字母（不要使用s c t  L l e）
       */
      double CalFunc(std::string &func, const std::string &arg, const std::string &arg_value);
//

      /*!
       * @brief 4th R-K Method
       * @details 4阶龙格库塔
       * @return 结果
       *    @retval double
       * @todo 未完成
       */
      static double RKMethod();  // 4th R-K Method

      /*!
       * @brief 公式法解三次方程
       * @details suppose equation to be ax^3+bx^2+cx+d=0？
       * @param a x^3 三次项系数
       * @param b x^2 二次项系数
       * @param c 一次项系数
       * @param d 常数
       * @return 解
       *    @retval double
       * @note 只有实根？
       */
      static double CubicEquRealRT(double a, double b, double c, double d);

      /*!
       * @brief 牛顿法解三次方程
       * @details suppose equation to be ax^3+bx^2+cx+d=0
       * @param x0 初始猜测值
       * @param err 误差额日日
       * @param a x^3 三次项系数
       * @param b x^2 二次项系数
       * @param c 一次项系数
       * @param d 常数
       * @return 误差范围内的解
       *    @retval double
       */
      static double NewtonCubicEquRT(double x0, double err, double a, double b, double c, double d);

  private:
      
      /*!
       * @brief 调度场算法转换为逆波兰式
       * @param a
       * @param b
       */
      static void Trans2Inverse(char *a, char *b);

      /*!
       * @brief 求解逆波兰式的值
       * @param expression
       * @return 逆波兰式的值
       *    @retval double
       */
      static double CompValue(char *expression);
  };


  
  /*!
   * @brief 图（数学意义）
   * 图的存储、处理和打印
   * 上面的描述中对术语“空图”的使用应该与空图或空图的数学定义区分开来
   * 严格地说，图论中的空图或空图是没有顶点和没有边的图
   * 但igraph所使用的“空图形”是指具有零个或多个顶点但没有边缘的图形
   */
  class Graph
  {
  public:
      static const int invalid = -9;

      /*!
       * @brief 最基本的构造函数，所有其他构造函数调用它来创建一个最小的图形对象
       */
      explicit Graph();

      /*!
       * @brief 由矩阵构造图
       * @param m 矩阵
       * @param directed 是否有向(true/false)
       * @attention 按行读取，不检查下标！请在从文件读入时注意控制溢出和方阵
       */
      explicit Graph(const std::vector< std::vector<double> > &m, bool directed);

      /*!
       * @brief 高维矩阵构造
       * @details 矩阵维数结构的向量 + 矩阵值
       * @param scale_vec 该向量大小决定网络的维数，每一个元素值分别代表该维度的格点数目
       * @param weight 权矩阵，按照(1,1,1,...,n) (2,1,1,...n) (3,1,1,...,n) (n,1,1,...,n) (n,2,1,...,n) ... 的顺序读取
       * 二维文本文件采用左上角坐标原点存储数据,如下
       *   1 2 3 4 ...
       * O------------>x
       * 1|
       * 2|
       * 3|
       * 4|
       * \/
       *  y
       */
      explicit Graph(std::vector<int> &scale_vec, std::vector<double> &weight);

      /*!
       * @brief 拷贝构造函数
       * @param G
       */
      Graph(const Graph& G);

      /*!
       * @brief 赋值运算符
       * @param G
       * @return 返回右值的引用
       */
      Graph& operator =(const Graph& G);

      /*!
       * @brief 析构运算
       */
      virtual ~Graph();
      
      /*!
       * @brief 是否有向图
       * @return 当前图对象是否有向(布尔值)
       *    @retval bool
       */
      inline bool IsDirected() const { return igraph_is_directed(&graf) != 0; }

      /*!
       * @brief 节点数
       * @return 当前图对象的节点数
       *    @retval int
       */
      inline int Size() const { return (int)igraph_vcount(&graf); }

      /*!
       * @brief 边数
       * @return 当前图对象的边数
       *    @retval int
       */
      inline int Edges() const { return (int)igraph_ecount(&graf); }
      
      /*!
       * @brief 改变图的向性
       * @param directed true有向 false无向
       * @return 改变后当前图的有向性
       *    @retval bool
       */
      bool ConvertDirectedness(bool directed);
      // 返回某一边的起点-终点对

      /*!
       * @brief 获得一条边的起点-终点对
       * @param edge_id 边的id
       * @return 起点和终点id组成的有序点对
       *    @retval std::pair<int,int>
       */
      std::pair<int,int> EdgeEndpoints(int edge_id) const;
      
      /*!
       * @brief 获得一对点间的有向（无向）边
       * @param vid_pair 起点和终点id组成的有序点对
       * @return 有向边的id，不存在边返回INVILID
       *    @retval int(INVILID)
       */
      int EdgeBetween( const std::pair<int, int> &vid_pair ) const;

      /*!
       * @brief 沿指定路径的边id序列
       * @details 图的向性会被考虑。可以用两种方式表示图的路径
       * @param vid_pairs
       * @param vid_along_path
       * @return 边id组成的向量
       *    @retval std::vector<int>
       */
      std::vector<int> EdgesAlongPath( const std::vector<std::pair<int, int>> &vid_pairs, const std::vector<int> &vid_along_path ) const;

      /*!
       * @brief 计算所有节点各自的出入度
       * @return 所有节点各自的出度和入度之和
       *    @retval std::vector<int>
       */
      std::vector<int> Degree() const;

      /*!
       * @brief 孤立节点统计
       * @return 图的孤立节点列表
       *    @retval std::vector<int>
       */
      std::vector<int> AloneVertices() const;

      /*!
       *
       * @return 是否强连通
       *    @retval bool
       */
      bool IsStrongConnected();

      /*!
       *
       * @return 是否弱联通
       *    @retval bool
       */
      bool IsWeakConnected();

      /*!
       *
       * @param vid_pair 起点-终点 id 组成的有序点对
       * @return 两点间的权值
       *    @retval double
       */
      inline double Weight( const std::pair<int, int> &vid_pair ) const{
          return MATRIX(wmat, vid_pair.first, vid_pair.second);
      }

      /*!
       *
       * @param vid_pairs
       * @return 一组点间的权值组成的向量
       *    @retval std::vector<double>
       */
      std::vector<double> Weight( const std::vector<std::pair<int, int>> &vid_pairs ) const{
          std::vector<double> res;
          for ( const auto &i : vid_pairs ){
            res.push_back(Weight(i));
          }
          return res;
      }
      
      /*!
       * @brief 向图中添加若干边
       * @param list 点对+权值组成的点对的向量
       * @return 添加是否成功
       *    @retval bool
       * @note 为什么不先实现添加一条边再添加多条边，因为igraph里添加多条边时，这样更快
       */
      bool AddEdges(std::vector<std::pair<std::pair<int, int> , double>> &list);

      /*!
       * @brief 向图中添加若干顶点
       * @param nv 想添加的点的数量
       * @return 添加是否成功
       *    @retval bool
       * @todo 同步权矩阵
       * @bug 未同步权矩阵
       */
      bool AddVertices(int nv = 1);

      /*!
       * @brief 向图中随机添加边
       * @param wlist
       * @return 是否成功
       *    @retval bool
       */
      bool RandAddEdges(std::vector<double> &wlist);

      /*!
       * @brief 删除图中边
       * @param vlist
       * @return 是否删除成功
       *    @retval bool
       */
      bool DelEdges(std::vector<std::pair<int, int>> &vlist);

      /*!
       * @brief 删除图中顶点
       * @param vid
       * @return 是否删除成功
       * @todo 同步权矩阵
       * @bug 未同步权矩阵
       */
      bool DelVertices(std::vector<int> &vid); // 注意！删除节点后序号将重新排列！！

      /*!
       * @brief 计算从起点到终点的无权最短路径
       * @param vid_pair 起点-终点
       * @return 无权最短路径值
       * @todo 用 igraph_get_shortest_path - 计算从/到一个顶点的最短路径 重写此函数
       */
      double UnweightShortPaths(const std::pair<int, int> &vid_pair) const;

      /*!
       * @brief 计算一组点间的无权最短路径
       * @param pair_vids 点对数组
       * @return 最短路径矩阵
       *    @retval std::vector<std::vector<double>>
       */
      std::vector<std::vector<double>> UnweightShortPaths(const std::pair<std::vector<int>, std::vector<int>> &pair_vids) const;

      /*!
       * @brief 计算所有两两点间的无权最短路径
       * @return 最短路径矩阵
       *    @retval std::vector<std::vector<double>>
       */
      std::vector<std::vector<double>> UnweightShortPaths() const;

      //!< @todo 还有一对多版本的 igraph_get_shortest_paths 也需要添加

      /*!
       * @brief 平均路径长度
       * @return 平均路径长度
       *    @retval double
       * @remark 似乎平均路径长度一直都是无权的
       */
      double AverPathLength() const;

      /*!
       * @brief 所有最短路径长度的分布/直方图
       * @return 直方图向量
       *    @retval std::vector<double>
       */
      std::vector<double> PathLengthDist() const;

      /*!
       * @brief 图的直径
       * @details 最长的测地线
       * @return 值
       *    @retval double
       * @todo 未实现
       * @bug 目前不可用
       */
      double Diameter() const; // 图的直径（最长测地线），未实现

      // 带权最短路径，两个均未实现
      /*!
       * @brief 带权最短路径
       * @param pair_vids 点向量组对
       * @return 最短路径矩阵
       *    @retval std::vector<std::vector<double>>
       * @todo 通过一个bool控制是否允许负权值？
       * @todo 未实现
       */
      std::vector<std::vector<double>> WeightShortPaths(const std::pair<std::vector<int>, std::vector<int>> &pair_vids) const;

      /*!
       * @brief 全局带权最短路径
       * @details 全部节点两两之间的带权最短路径
       * @param pair_vids 点向量组对
       * @return 最短路径矩阵
       *    @retval std::vector<std::vector<double>>
       * @todo 通过一个bool控制是否允许负权值？
       * @todo 未实现
       */
      std::vector<std::vector<double>> WeightShortPaths() const;

      /*!
       * @brief n近临节点
       * @details 和一个节点数组（里面是节点vid）里的节点n近临的节点组成的矩阵
       * The neighborhood of a given order of a vertex includes all vertices which are closer to the vertex than the order.
       * Ie. order 0 is always the vertex itself, order 1 is the vertex plus its immediate neighbors,
       * order 2 is order 1 plus the immediate neighbors of the vertices in order 1, etc.
       * This function calculates the vertices within the neighborhood of the specified vertices.
       * 
       * @param vids
       * @param order
       * @return 节点矩阵
       *    @retval std::vector<std::vector<int>>
       * @todo 未实现
       */
      std::vector<std::vector<int>> Neighborhood(std::vector<int> vids, int order) const;

      /*!
       * @brief 图的铰接点
       * @details 如果顶点移除增加了图中连接组件的数量，则顶点是铰接点
       * @return 铰接点的枪毙名单
       *    @retval int
       */
      std::vector<int> ArticulationPoints() const;


      // igraph_closeness

      // igraph_betweenness

      // igraph_edge_betweenness

      // igraph_strength - 顶点的强度，换句话说加权顶点度

      /*!
       * @brief 全局最大度
       * @details Calculate the maximum degree of the graph
       * @return 一个最大的度值
       *    @retval int
       */
      int MaxDegree() const;

      /*!
       * @brief 局部最大度
       * @details Calculate the maximum degree of a set of vertices in the graph
       * @return 一个最大的度
       *    @retval int
       */
      int MaxDegree(std::vector<int> vids) const;

      /*!
       * @brief 最大度点集
       * @return 具有最大度的点的列表
       *    @retval std::vector<int>
       */
      std::vector<int> MaxDegreePoints() const;

      // 中心度
      // igraph_centralization_degree
      // igraph_centralization_betweenness

      /*!
       * @brief 全局聚集系数
       * @return 全局聚集系数
       *    @retval double
       */
      double ClusterCoeff() const;
      // 全局的聚集系数
      // 局部的聚集系数？
      // 加权的聚集系数？ igraph_transitivity_barrat

      // 是否是简单图，删除循环等等

      /*!
       * @brief 最小分割集
       * Find all minimum size separating vertex sets
       * This function lists all separator vertex sets of minimum size.
       * A vertex set is a separator if its removal disconnects the graph.
       *
       * @return 点集矩阵
       *    @retval std::vector<std::vector<int>>
       */
      std::vector<std::vector<int>> MiniumSizeSeps() const;

      // K-Cores ?

  protected:
      //!< @todo 利用igraph中的一些生成经典图方法，实现某一具体的网络。只能被子类调用，不开放对外接口

  private:
      igraph_t graf;        //!< igraph 图对象（结构）
      igraph_matrix_t wmat; //!< igraph 权矩阵

      /*!
       * @brief 更新权值
       * @details 更新指定节点间的权值。有向图表示从 i 到 j，无向图表示 i 和 j 之间
       * @param i 从 i
       * @param j 到 j
       * @param new_w 新的权值
       * @attention 更新权值时注意考虑图的向性
       */
      void UpdateWMat(int i, int j, double new_w){
          igraph_matrix_set(&wmat, i, j, new_w);
          if(!igraph_is_directed(&graf)) {
              igraph_matrix_set(&wmat, i, j, new_w);
          }
      }

      /*!
       * @brief 将 STL 向量转化为 igraph 向量
       * @tparam T
       * @param std_vec 标准库向量
       * @param igraph_vec 未初始化的 igraph_vector_t & 变量
       * @warning 需要提供未初始化的igraph_vector_t指针，引用？
       */
      template <typename T>
      void Vec2igVec(const std::vector<T> &std_vec, igraph_vector_t &igraph_vec) const;

      /*!
       * @brief 将 igraph 向量转化为 STL 向量
       * @tparam T
       * @param igraph_vec
       * @return 转化结果（STL向量）
       *    @retval std::vector<T>
       */
      template <typename T>
      std::vector<T> IgVec2intVec(const igraph_vector_t &igraph_vec) const;

      /*!
       * @brief 将 pair 组成的 vector 拍扁成 igraph 向量
       * @tparam T
       * @param src 作为来源的 std::vector<std::pair<T, T>> &src
       * @param dst 未初始化的 igraph_vector_t & 变量
       * @warning 需要提供未初始化的igraph_vector_t指针，引用？
       */
      template <typename T>
      void FlatVecPairs2igVec(const std::vector<std::pair<T, T>> &src, igraph_vector_t &dst) const;

      /*!
       * @brief 将一个 pair 组成的数组切分为两个 igraph 数组
       * @tparam T
       * @param src std::vector<std::pair<T, T>> &
       * @param dst1 igraph_vector_t &
       * @param dst2 igraph_vector_t &
       * @warning 需要提供未初始化的igraph_vector_t指针，引用？
       */
      template <typename T>
      void DivVecPairs2igVecs(const std::vector<std::pair<T, T>> &src,
                              igraph_vector_t &dst1, igraph_vector_t &dst2) const;

      /*!
       * @brief 将两个 vector 数组组成的 pair 切分为两个 igraph 数组
       * @tparam T
       * @param std::pair<std::vector<T>, std::vector<U>> &
       * @param dst1 igraph_vector_t &
       * @param dst2 igraph_vector_t &
       * @warning 需要提供未初始化的igraph_vector_t指针，引用？
       */
      template <typename T, typename U>
      void DivPairVecs2igVecs(const std::pair<std::vector<T>, std::vector<U>> &src,
                              igraph_vector_t &dst1, igraph_vector_t &dst2 ) const;

      /*!
       * @brief 将 igraph 的矩阵转换为 vec * vec 矩阵
       * @tparam T
       * @param src igraph_matrix_t &
       * @return vec * vec 矩阵
       *    @retval std::vector<std::vector<T>>
       */
      template<typename T>
      std::vector<std::vector<T>>IgMatrix2VecVecs(const igraph_matrix_t &src) const;

      
      template <typename T, typename U>
      std::vector<T> ExtVecPairs1stVec(const std::vector<std::pair<T,U>> &src) const;
      template <typename T, typename U>
      std::vector<U> ExtVecPairs2ndVec(const std::vector<std::pair<T,U>> &src) const;
  };

}   // namespace ntoolkit

#endif //NTOOLKIT_VLIB_H
