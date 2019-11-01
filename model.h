//
// Created by verniy on 18-10-29.
//

#ifndef NTOOLKIT_MODEL_H
#define NTOOLKIT_MODEL_H

#define LINUX
//#define WIN

#include <vector>

#ifdef LINUX
#include <tbb/parallel_for.h>
#include <tbb/task.h>
#elif WIN
#include "tbb/parallel_for.h"
#include "tbb/task.h"
#endif

#include "neuron.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_traits.hpp>  // 顶点迭代器类型

namespace std {
    template <typename T>
    std::istream& operator>>(std::istream& in, std::pair<T,T>& p) {
        in >> p.first >> p.second;
        return in;
    }
    template <typename T, typename U>
    std::istream& operator>>(std::istream& in, std::pair<std::pair<T,T>, U>& p) {
        in >> p.first.first >> p.first.second >> p.second;
        return in;
    }
}

namespace boost {
    enum vertex_neuron_obj_t { vertex_neuron_obj = 111 };
    enum vertex_coordinate_t { vertex_coordinate = 112 };
    BOOST_INSTALL_PROPERTY(vertex, neuron_obj);
    BOOST_INSTALL_PROPERTY(vertex, coordinate);

        /*
    提供给boost架构，用boost来实现下面所有功能
    pair<int,int> EdgeEndpoints(int edge_id) const;    // 返回某一边的起点-终点对
    int EdgeBetween( const pair<int, int> &vid_pair ) const;
    // 起/终点间是否存在边（考虑图的向性），返回这个边的id，不存在返回INVILID
    vector<int> EdgesAlongPath( const vector<pair<int, int>> &vid_pairs, const vector<int> &vid_along_path ) const;
    // 沿指定路径的边id序列（考虑图的向性）

    vector<int> Degree() const;         // 返回所有节点的出度和入度之和
    vector<int> AloneVertices() const;  // 返回图的孤立节点列表
    bool IsStrongConnected();           // 强连通？
    bool IsWeakConnected();             // 弱连通？
    double Weight( const pair<int, int> &vid_pair ) const { return MATRIX(wmat, vid_pair.first, vid_pair.second); }        // 返回两点间权值
    vector<double> Weight( const vector<pair<int, int>> &vid_pairs ) const;                        // 返回一组点间的权值

    // 未添加权矩阵的更新！！
    bool AddEdges(vector<pair<pair<int, int> , double>> &list);
    bool AddVertices(int nv = 1);       // 默认添加一个点，修改矩阵的功能未添加！！！
    bool RandAddEdges(vector<double> &wlist);   // 随机添加连接
    bool DelEdges(vector<pair<int, int>> &vlist);
    bool DelVertices(vector<int> &vid); // 注意！删除节点后序号将重新排列！！
    // 修改矩阵的功能未添加！！！

    double UnweightShortPaths(const pair<int, int> &vid_pair) const;
    // 请用 igraph_get_shortest_path - 计算从/到一个顶点的最短路径 重写此函数
    vector<vector<double>> UnweightShortPaths(const pair<vector<int>, vector<int>> &pair_vids) const;
    vector<vector<double>> UnweightShortPaths() const;   // 返回所有点两两配对组合的最短路径，是为矩阵
    // 还有一对多版本的 igraph_get_shortest_paths 也需要添加

    double AverPathLength() const; // 似乎平均路径长度一直都是无权的
    vector<double> PathLengthDist() const; // 所有最短路径长度的分布/直方图
    double Diameter() const; // 图的直径（最长测地线）

    vector<vector<double>> WeightShortPaths(const pair<vector<int>, vector<int>> &pair_vids) const; // 加权最短路径，通过一个bool控制是否允许负权值？
    vector<vector<double>> WeightShortPaths() const;
    // 未实现

    vector<vector<int>> Neighborhood(vector<int> vids, int order) const;
     The neighborhood of a given order of a vertex includes all vertices which are closer to the vertex than the order.
     * Ie. order 0 is always the vertex itself, order 1 is the vertex plus its immediate neighbors, order 2 is order 1 plus
     *
     * the immediate neighbors of the vertices in order 1, etc.
     * This function calculates the vertices within the neighborhood of the specified vertices.
     *

    vector<int> ArticulationPoints() const;
    // 如果顶点的移除增加了图中连接组件的数量，则顶点是铰接点

    // igraph_closeness

    // igraph_betweenness

    // igraph_edge_betweenness

    // igraph_strength - 顶点的强度，换句话说加权顶点度

    int MaxDegree() const;
    // Calculate the maximum degree in a graph
    int MaxDegree(vector<int> vids) const;
    // in a set of vertices

    vector<int> MaxDegreePoints() const;
    // 返回具有最大度的点的列表

    // igraph_centralization_degree
    // igraph_centralization_betweenness

    double ClusterCoeff() const;
    // 全局的聚集系数
    // 局部的聚集系数？
    // 加权的聚集系数？ igraph_transitivity_barrat

    // 是否是简单图，删除循环等等

    vector<vector<int>> MiniumSizeSeps() const;
    // Find all minimum size separating vertex sets
    // This function lists all separator vertex sets of minimum size.
    // A vertex set is a separator if its removal disconnects the graph.

    // K-Cores ?

    protected:
    // 利用igraph中的一些生成经典图的方法，但只能被子类调用（用来实现某一个具体的网络），不开放对外接口

    private:


    virtual AddVetex();
    */

    // boost域里加入我自己的模板类，模板类声明定义写在一起

    template <typename neuron_t>
    class Neuron_net
    {
    public:
        Neuron_net() = default;
        explicit Neuron_net(const std::string &mat_file_path, bool is_sparse);
        unsigned int Size() const { return boost::num_vertices(G); }
        unsigned int EdgeNum() const { return boost::num_edges(G); }
        std::string Name (unsigned int id) const {
            vertex_it_t it, it_end;
            //index_map_t vid = boost::get(boost::vertex_index, G);
            //name_map_t vname = boost::get(boost::vertex_name, G);
            for(std::tie(it, it_end) = boost::vertices(G); it != it_end; ++it){
                if(boost::get(boost::vertex_index, G, *it) == id)
                    return boost::get(boost::vertex_name, G, *it);
            }
        }

    private:
        typedef boost::adjacency_matrix<
                boost::vecS,
                boost::vecS,
                boost::undirectedS,
                boost::property< boost::vertex_index_t, int,
                        boost::property<boost::vertex_name_t, std::string,
                        boost::property<boost::vertex_coordinate_t, std::tuple<double, double, double>,
                        boost::property<boost::vertex_neuron_obj_t, neuron_t > > > >,
                boost::property<boost::edge_weight_t, double>
        > graph_t;
        typedef typename boost::property_map<graph_t, boost::vertex_index_t>::type index_map_t;
        typedef typename boost::property_map<graph_t, boost::vertex_name_t>::type name_map_t;
        typedef typename boost::property_map<graph_t, boost::vertex_coordinate_t>::type coordinate_map_t;
        typedef typename boost::property_map<graph_t, boost::vertex_neuron_obj_t>::type neuron_obj_map_t;

        //typedef typename boost::graph_traits<graph_t>::adjacency_iterator graph_it_t;
        typedef typename boost::graph_traits<graph_t>::vertex_iterator vertex_it_t;
        typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_dp_t;

        graph_t G;
    };

    template <typename neuron_t>
    Neuron_net<neuron_t>::Neuron_net(const std::string &mat_file_in_path, bool is_sparse)
    {
        if(is_sparse) {
            std::ifstream file_in(mat_file_in_path);

            typedef typename boost::graph_traits<graph_t>::vertices_size_type size_type;
            size_type n_vertices;   // 文件的第一个元素标识节点数量
            file_in >> n_vertices; // read in number of vertices
            // 如何有效处理权值属性
            typedef std::pair<std::pair<size_type, size_type>, boost::property<boost::edge_weight_t, double> > edge_t;
            std::istream_iterator<edge_t> file_in_it(file_in), eof;
            std::vector<edge_t> file_in_vec(file_in_it, eof);
            // 读入
            std::vector< std::pair<size_type, size_type> > edges_vec; edges_vec.clear();
            std::vector< boost::property<boost::edge_weight_t, double> > weights_vec; weights_vec.clear();
            for(const auto &it : file_in_vec){
                edges_vec.push_back(it.first);
                weights_vec.push_back(it.second);
            }
            // 分成两个数组，不知道有没有什么更好的办法

            G(edges_vec.begin(), edges_vec.end(), n_vertices);
            // 是否有更好的处理办法

        } else {

        }
    }
}

namespace ntoolkit{

    class Neuron_net_base
    {
    public:
        Neuron_net_base();
        virtual ~Neuron_net_base();
        virtual std::vector<double> Refresh();
        virtual std::vector<double> Refresh(const std::vector<double> &ext);

        ntoolkit::math::Graph* Network() const { return net; }

        bool GenNetReport(std::string filename);
        // 生成网络报告到文件
    protected:
        ntoolkit::math::Graph* net;
        std::vector<ntoolkit::Neuron*> neuron;
    private:
        void Clear();
    };


    template < typename Neuron_T, typename Net_T > class Neuron_net : public Neuron_net_base {
    public:
        Neuron_net(const std::vector<double> &initial, const std::vector<std::vector<double>> &w, std::initializer_list<double> parameter);
    };

    template < typename Neuron_T, typename Net_T >
    Neuron_net<Neuron_T, Net_T>::Neuron_net(const std::vector<double> &initial, const std::vector<std::vector<double>> &w, std::initializer_list<double> parameter)
    {
        if (initial.empty()) throw "初始向量不能为空";
        if (w.size() != initial.size()) throw "初始向量和权矩阵维度不符";
        for (const double& i : initial){
            ntoolkit::Neuron* temp;
            temp = new Neuron_T(i,parameter);
            neuron.push_back(temp);
        }
        net = new Net_T(w);
    }

/*
class AiharaBrainNet : public Neuron_net
{
public:
    AiharaBrainNet();
    void InitialNeuron(const vector<double> &paras_name_list, const vector<double> &initial) override;
    void InitialNet(const double *weight) override;
};


class HHBrainNet : public Neuron_net
{
public:
    HHBrainNet();
    void InitialNeuron(const vector<double> &paras_name_list, const vector<double> &initial) override;
    void InitialNet(const double *weight) override;
};


class FHNBrainNet : public Neuron_net
{
public:
    FHNBrainNet();
    void InitialNeuron(const vector<double> &paras_name_list, const vector<double> &initial) override;
    void InitialNet(const double *weight) override;
};

class FHNBrainNet2 : public Neuron_net
{
public:
    FHNBrainNet2();
    void InitialNeuron(const vector<double> &paras_name_list, const vector<double> &initial) override;
    void InitialNet(const double *weight) override;
};
*/
}

#endif //NTOOLKIT_MODEL_H
