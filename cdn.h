//
//  cdn.h
//  CDN复赛
//
//  Created by Daniel on 2017/4/9.
//  Copyright © 2017年 Daniel. All rights reserved.
//

#ifndef cdn_h
#define cdn_h

#include <stdio.h>
#include <stdlib.h>

#include "deploy.h"
#include "lib_io.h"
#include "lib_time.h"

#include <bitset>
#include <time.h>
#include <unordered_set>
#include <vector>
#include <queue>
#include <iostream>

#include <algorithm>
#include <sys/timeb.h>
#include <cmath>
#include <fstream>
#include <map>
#include <cstdlib>

#include <list>
#include <string>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <unordered_map>

#include <utility>
#include <set>


#define TIME_OUT 88 // 85s



#define MAX_NUMBER 1215752191

#define MAX(a,b) (((a)>(b))?a:b)
#define MIN(a,b) (((a)<(b))?a:b)

#define MaxNetNode          10000                // 网络节点上限
#define MaxConsumNode       10000                 // 消费节点上限

#define IB_BOTTLENECK_ORIG 0
#define IBTEST 0
#define IB_MIN_MARGINALS_DEBUG 0
#define IB_MIN_MARGINALS_TEST 0
#define IBSTATS 0
#define IBDEBUG(X) fprintf(stdout, "\n"); fflush(stdout)
#define IB_ALTERNATE_SMART 1
#define IB_HYBRID_ADOPTION 1
#define IB_EXCESSES 1
#define IB_ALLOC_INIT_LEVELS 4096
#define IB_ADOPTION_PR 0
#define IB_DEBUG_INIT 0
#define IB_PREVPTR_EXCESS(x) (ptrs[(((x)-nodes)<<1) + 1])
#define IB_NEXTPTR_EXCESS(x) (ptrs[((x)-nodes)<<1])
#define IB_PREVPTR_3PASS(x) ((x)->firstSon)

#define REMOVE_SIBLING(x, tmp) \
{ (tmp) = (x)->parent->head->firstSon; \
if ((tmp) == (x)) { \
(x)->parent->head->firstSon = (x)->nextPtr; \
} else { \
for (; (tmp)->nextPtr != (x); (tmp) = (tmp)->nextPtr); \
(tmp)->nextPtr = (x)->nextPtr; \
} }

#define ADD_SIBLING(x, parentNode) \
{ (x)->nextPtr = (parentNode)->firstSon; \
(parentNode)->firstSon = (x); \
}

using namespace std;


typedef struct bfs{
    friend bool operator < (bfs n1, bfs n2){
        return n1.pathcost < n2.pathcost;
    }
    vector<int> node;
    int pathcost=0;
    
} bfsnode;


class ServerNode{
public:
    int netnode;    // 网络节点id
    int serverLevel;// 服务器档次
    int flow;
    ServerNode(){}
    ServerNode(int a,int b) {
        netnode = a;
        serverLevel = b;
        flow = 0;
    }
    ServerNode(int a,int b,int c) {
        netnode = a;
        serverLevel = b;
        flow = c;
    }
    
};

struct NetNode
{
    friend bool operator < (NetNode n1, NetNode n2){
        return n1.flow_ <= n2.flow_;
    }
    
    NetNode() : cost(0), flow_(0),  degree_(0) {}
    
    int cost = 0; // 部署成本
    int flow_ = 0; // 网络节点总流量
    int degree_ = 0; // 度
    int maxServerLevel = 0; // 按输能力，合适的档次
    float unit_cost = 0; // 带宽单位成本
    
    void Print() {
        cout<<"f:"<<flow_<<" d:"<<degree_<<endl;
    }
};
struct ConsumNode{
    int netNode;
    int need;
    
    ConsumNode() : netNode(0), need(0) {}
    
};
struct ServerLevel{
    int levelid;
    int outPut;
    int cost;
    
    ServerLevel() : levelid(0), outPut(0),  cost(0) {}
};

class NetLink {
public:
    int startPoint;     // 链路起点
    int endPoint;       // 链路终点
    int maxWidth;       // 总带宽
    int remainWidth;    // 剩余的宽带
    int cost;           // 单位带宽租金
    NetLink* next;
    
    NetLink() {}
    NetLink(int s,int e,int w,int c) {
        this->startPoint = s;
        this->endPoint = e;
        this->maxWidth = w;
        this->cost = c;
        this -> next = NULL;
    }
};

struct RouteDetail
{
    vector<ServerNode> server_set;   // 服务器节点集合 带宽单位成本 从大到小排序 
    vector<ServerNode> pass_set;    // 路由节点集合  带宽单位成本 从小到大排序
    vector<int> up_set;    // 升档集 服务器在路由路径中
    unordered_map<int, unordered_set<int>> route_set;    // 服务器的路由节点集
    unordered_map<int, unordered_set<int>> service_set;    // 服务器服务的消费节点集 
    /* 
    type: 1 server_set 删点&降挡
    type: 2 pass_set   松弛加点
    type: 3 up_set 升档
    type: 4 server_set route_set 服务器路由节点集合 交换
    */
    int type;     
    
};

typedef NetLink* Graph;


class TopoInfo {
public:
    vector<Graph> graph;                    // 邻接表 存图
    vector<NetLink> g_edges_maxflow;        // 最大流需要用到的边
    vector<int> g_net2consume;             // 网络节点对消费者节点
    int NetNodeCount;                       // 总网络节点个数
    int LinkCount;                          // 总链路条数
    int ConsumNodeCount;                    // 总消费节点个数
    int ServerLevelCount;                   
    vector<ServerLevel> serverLevel;        // 服务器档次
    vector<ConsumNode> consumNode;
    vector<NetNode> net_nodes;               // 网络顶点表
    bitset<MaxNetNode> key_bit;
    int all_demend; //总需求
    int width_unit;     // 带宽单位成本 
};



typedef std::pair<int,int>P;

class Flow {
public:
    int cur_node;
    int need;
    string road;
    vector<int> roads;
};


class Edge{
public:
    int start,end,flow,rent; // s：边尾, t: 边头 c: 流量 r: 租金
    Edge(int a,int b,int c):start(a),end(b),flow(c) {}
    Edge(int a,int b,int c, int d):start(a),end(b),flow(c),rent(d) {}
    Edge(){}
};


class McmfEdge;

class McmfNode {
public:
    long excess; // 节点流量流出边
    int price;  // 到汇点的距离
    McmfEdge *first; // 第一条邻接边
    McmfEdge *current; // 当前指向临接边
    McmfEdge *suspended; // 挂起边
    McmfNode *q_next; // 队列里的下一个节点
    McmfNode *b_next; // 桶里的下一个节点
    McmfNode *b_prev; // 桶里的上一个节点
    int rank;    // 桶数量
    int inp;     // 临时边的数量
    int cici;
    
    McmfNode() {}
    ~McmfNode() {}
};

class McmfEdge {
public:
    int rez_capacity; // 剩余容量
    int cost; // 边的成本
    McmfNode *head;	 // 边的头节点
    McmfNode *tail;	// 边的尾节点
    McmfEdge *sister; // 相反的边
    int sisi;
    
    McmfEdge() {}
    ~McmfEdge() {}
};

class McmfBucket {
public:
    McmfNode *p_first; 	// 第一个正出量的节点或者桶的第一个节点
    McmfEdge *ee;
    
    McmfBucket( McmfNode *p_first) : p_first(p_first) {}
    McmfBucket() {}
    ~McmfBucket() {}
};


class Mcmf {
public:
    vector<Edge> roads; // 路由流量记录
    int cost_maxflow; // 最小费用流
    unordered_map<int, int> server_supply; // 服务器流量提供量
    
    bool is_error; // 全局错误状态
    
    int n; // 节点数量
    int m; // 边的数量
    
    int *cap; // 边初始流量数组
    McmfNode *nodes; // 节点数组
    McmfNode *sentinel_node; // 最后节点的下一个
    McmfNode *excq_first; // 队列的队头
    McmfNode *excq_last; // 队列的队尾
    McmfEdge *arcs; // 边的数组
    McmfEdge *sentinel_arc; //最后一条边的下一条边
    
    McmfBucket *buckets; // 桶的数组
    McmfBucket *l_bucket; // 最后一个桶
    int linf; // 桶数量+1
    int time_for_price_in;
    
    int epsilon; // 最优化约束
    int dn;     // cost乘数
    int price_min; // 下限
    int mmc;  // 最大成本
    double f_scale; // 比例因子
    double cut_off_factor;
    double cut_on; //返回挂起边的约束
    double cut_off; // 挂起边的约束
    long total_excess; // 出量总数
    
    int flag_price;
    int flag_updt;
    
    McmfEdge d_arc; // 傀儡边
    McmfNode *dnode;
    
    int n_rel; // 重标签数量
    int n_ref; // 当前改善的数量
    int n_src;  // 当前超量的节点数量
    int n_bad_pricein;
    int n_bad_relabel;
    
    int node_min; // 节点最小序号
    int node_max; // 节点最大序号
    int *arc_first; // 节点度数组
    
    int *arc_tail; // 边尾数组
    int pos_current; // 记录当前指向节点
    McmfEdge *arc_current;  // 记录当前指向边
    McmfEdge *arc_new; // 新的边
    McmfEdge *arc_tmp; // 临时边
    int max_cost; // 最大权值
    int total_p; // 提供总流量
    int total_n; // 总需求
    
    McmfNode *i_node; // 添加边使用
    McmfNode *j_node; // 添加边使用
    
    int qiqi;
    
public:
    Mcmf(int num_nodes, int num_arcs) {
        is_error = false;
        cost_maxflow = -1.0;
        n = num_nodes;
        m = num_arcs;
        
        node_max = n - 1;
        node_min = 0;
        flag_price = 0;
        flag_updt = 0;
        n_bad_pricein = 0;
        n_bad_relabel = 0;
        
        allocate_arrays();
    }
    ~Mcmf() {}
    
    void allocate_arrays();
    void deallocate_arrays();
    void set_arc(int tail_id, int head_id, int capacity,int cost);
    void set_supply_demand_of_node( int head_id, int excess);
    void pre_processing();
    void algo_initialize();
    void up_node_scan( McmfNode *i);
    void price_update();
    int relabel( McmfNode *i);
    void discharge( McmfNode *i);
    int price_in();
    void refine();
    int price_refine();
    void price_out();
    int update_epsilon();
    void build_result(double *objective_cost);
    void algo( double *objective_cost);
    int exec();
    
    // 共享utils
    void increase_flow( McmfNode *i, McmfNode *j, McmfEdge *a, long df);
    // 超量队列的utils
    void reset_excess_q();
    void insert_to_excess_q( McmfNode *i) ;
    void insert_to_bucket( McmfNode *i, McmfBucket *b) ;
    void remove_from_bucket( McmfNode *i, McmfBucket *b);
    void update_cut_off();
    void exchange( McmfEdge *a, McmfEdge *b);
};


/**
 图的初始化
 
 @param topo 图文件
 @param line_num 行数
 @param topoInfo 图的基本信息
 @return 0
 */
int initial(char * topo[MAX_EDGE_NUM], int line_num, TopoInfo& topoInfo);




/**
 时间检查
 
 @return 超时返回ture
 */
bool TimeOut();



/**
 输出用最大流获取的结果
 
 @param result 存放输出字符串
 @param g_servers 服务器
 @param g_best_roads 使用到的链路集合
 @param NetNodeCount 网络节点个数
 @param g_net2consume 网络节点id转消费节点id
 */
void writeMaxFlowPrograme(string& result, vector<ServerNode>& g_servers, vector<Edge>& g_best_roads, int NetNodeCount, vector<int>& g_net2consume);


/**
 检查适合的服务器档次

 @param output 服务器输出量
 @return 服务器档次
 */
ServerLevel checkServerLevel(int output);

class MfGraph
{
public:
    vector<int> fail_consumer;

    MfGraph():prNodeBuckets(orphan3PassBuckets)
    {
        arcIter = NULL;
        incList = NULL;
        incLen = incIteration = 0;
        numNodes = 0;
        uniqOrphansS = uniqOrphansT = 0;
        augTimestamp = 0;
        verbose = IBTEST;
        arcs = arcEnd = NULL;
        nodes = nodeEnd = NULL;
        topLevelS = topLevelT = 0;
        flow = 0;
        memArcs = NULL;
        tmpArcs = NULL;
        tmpEdges = tmpEdgeLast = NULL;
        ptrs = NULL;
    }


    ~MfGraph()
    {
        delete []nodes;
        delete []memArcs;
        orphanBuckets.free();
        orphan3PassBuckets.free();
        excessBuckets.free();
    }

    void setVerbose(bool a_verbose) {
        verbose = a_verbose;
    }

    void initSize(int numNodes, int numEdges);
    void addEdge(int nodeIndexFrom, int nodeIndexTo, int capacity, int reverseCapacity);
    void addNode(int nodeIndex, int capFromSource, int capToSink);
    void incEdge(int nodeIndexFrom, int nodeIndexTo, int capacity, int reverseCapacity);
    void incNode(int nodeIndex, int deltaCapFromSource, int deltaCapToSink);
    bool incShouldResetTrees();
    struct Arc;
    void incArc(Arc *a, int deltaCap);
    void initGraph();
    int computeMaxFlow();
    int computeMaxFlow(bool allowIncrements);
    void resetTrees();
    void computeMinMarginals();
    void pushRelabel();
    void analysis_fail();

    inline int getFlow() {
        return flow;
    }
    inline int getNumNodes() {
        return (int)(nodeEnd-nodes);
    }
    inline int getNumArcs() {
        return (int)(arcEnd-arcs);
    }
    int isNodeOnSrcSide(int nodeIndex, int freeNodeValue = 0);

    struct Node;

    struct Arc
    {
        Node* head;
        Arc*    rev;
        int isRevResidual :1;
        int rCap :31;
    };

    struct Node
    {
        int lastAugTimestamp:30;
        int isParentCurr:1;
        int isIncremental:1;
        Arc *firstArc;
        Arc *parent;
        Node *firstSon;
        Node *nextPtr;
        int label;
        int excess;
    };

private:
    Arc *arcIter;
    void augment(Arc *bridge);
    template<bool sTree> int augmentPath(Node *x, int push);
    template<bool sTree> int augmentExcess(Node *x, int push);
    template<bool sTree> void augmentExcesses();
    template<bool sTree> void augmentDischarge(Node *x);
    template<bool sTree> void augmentExcessesDischarge();
    template<bool sTree> void augmentIncrements();
    template <bool sTree> void adoption(int fromLevel, bool toTop);
    template <bool sTree> void adoption3Pass(int minBucket);
    template <bool dirS> void growth();

    int computeMaxFlow(bool trackChanges, bool initialDirS);
    void resetTrees(int newTopLevelS, int newTopLevelT);

    template<bool sTree> void pushRelabelDischarge(Node *x);
    template<bool sTree> void pushRelabelGlobalUpdate();
    template<bool sTree> void pushRelabelDir();
    void pushRelabelShelve(int fromLevel);

    class ActiveList
    {
    public:
        inline ActiveList() {
            list = NULL;
            len = 0;
        }
        inline void init(Node **mem) {
            list = mem;
            len = 0;
        }
        inline void clear() {
            len = 0;
        }
        inline void add(Node* x) {
            list[len] = x;
            len++;
        }
        inline Node* pop() {
            len--;
            return list[len];
        }
        inline Node** getEnd() {
            return list+len;
        }
        inline static void swapLists(ActiveList *a, ActiveList *b) {
            ActiveList tmp = (*a);
            (*a) = (*b);
            (*b) = tmp;
        }
        Node **list;
        int len;
    };

    class BucketsOneSided
    {
    public:
        inline BucketsOneSided() {
            buckets = NULL;
            maxBucket = 0;
            nodes = NULL;
            allocLevels = 0;
        }
        inline void init(Node *a_nodes, int numNodes) {
            nodes = a_nodes;
            allocLevels = numNodes/8;
            if (allocLevels < IB_ALLOC_INIT_LEVELS) {
                if (numNodes < IB_ALLOC_INIT_LEVELS) allocLevels = numNodes;
                else allocLevels = IB_ALLOC_INIT_LEVELS;
            }
            buckets = new Node*[allocLevels+1];
            memset(buckets, 0, sizeof(Node*)*(allocLevels+1));
            maxBucket = 0;
        }
        inline void allocate(int numLevels) {
            if (numLevels > allocLevels) {
                allocLevels <<= 1;
                Node **alloc = new Node*[allocLevels+1];
                memset(alloc, 0, sizeof(Node*)*(allocLevels+1));
                delete []buckets;
                buckets = alloc;
            }
        }
        inline void free() {
            delete []buckets;
            buckets = NULL;
        }
        template <bool sTree> inline void add(Node* x) {
            int bucket = (sTree ? (x->label) : (-x->label));
            x->nextPtr = buckets[bucket];
            buckets[bucket] = x;
            if (bucket > maxBucket) maxBucket = bucket;
        }
        inline Node* popFront(int bucket) {
            Node *x;
            if ((x = buckets[bucket]) == NULL) return NULL;
            buckets[bucket] = x->nextPtr;
            return x;
        }

        Node **buckets;
        int maxBucket;
        Node *nodes;
        int allocLevels;
    };
    class Buckets3Pass
    {
    public:
        inline Buckets3Pass() {
            buckets = NULL;
            nodes = NULL;
            maxBucket = allocLevels = -1;
        }
        inline void init(Node *a_nodes, int numNodes) {
            nodes = a_nodes;
            allocLevels = numNodes/8;
            if (allocLevels < IB_ALLOC_INIT_LEVELS) {
                if (numNodes < IB_ALLOC_INIT_LEVELS) allocLevels = numNodes;
                else allocLevels = IB_ALLOC_INIT_LEVELS;
            }
            buckets = new Node*[allocLevels+1];
            memset(buckets, 0, sizeof(Node*)*(allocLevels+1));
            maxBucket = 0;
        }
        inline void allocate(int numLevels) {
            if (numLevels > allocLevels) {
                allocLevels <<= 1;
                Node **alloc = new Node*[allocLevels+1];
                memset(alloc, 0, sizeof(Node*)*(allocLevels+1));
                delete []buckets;
                buckets = alloc;
            }
        }
        inline void free() {
            delete []buckets;
            buckets = NULL;
        }
        template <bool sTree> inline void add(Node* x) {
            int bucket = (sTree ? (x->label) : (-x->label));
            if ((x->nextPtr = buckets[bucket]) != NULL) IB_PREVPTR_3PASS(x->nextPtr) = x;
            buckets[bucket] = x;
            if (bucket > maxBucket) maxBucket = bucket;
        }
        inline Node* popFront(int bucket) {
            Node *x = buckets[bucket];
            if (x == NULL) return NULL;
            buckets[bucket] = x->nextPtr;
            IB_PREVPTR_3PASS(x) = NULL;
            return x;
        }
        template <bool sTree> inline void remove(Node *x) {
            int bucket = (sTree ? (x->label) : (-x->label));
            if (buckets[bucket] == x) {
                buckets[bucket] = x->nextPtr;
            } else {
                IB_PREVPTR_3PASS(x)->nextPtr = x->nextPtr;
                if (x->nextPtr != NULL) IB_PREVPTR_3PASS(x->nextPtr) = IB_PREVPTR_3PASS(x);
            }
            IB_PREVPTR_3PASS(x) = NULL;
        }
        inline bool isEmpty(int bucket) {
            return buckets[bucket] == NULL;
        }

        Node **buckets;
        int maxBucket;
        Node *nodes;
        int allocLevels;
    };

    class ExcessBuckets
    {
    public:
        inline ExcessBuckets() {
            buckets = ptrs = NULL;
            nodes = NULL;
            allocLevels = maxBucket = minBucket = -1;
        }
        inline void init(Node *a_nodes, Node **a_ptrs, int numNodes) {
            nodes = a_nodes;
            allocLevels = numNodes/8;
            if (allocLevels < IB_ALLOC_INIT_LEVELS) {
                if (numNodes < IB_ALLOC_INIT_LEVELS) allocLevels = numNodes;
                else allocLevels = IB_ALLOC_INIT_LEVELS;
            }
            buckets = new Node*[allocLevels+1];
            memset(buckets, 0, sizeof(Node*)*(allocLevels+1));
            ptrs = a_ptrs;
            reset();
        }
        inline void allocate(int numLevels) {
            if (numLevels > allocLevels) {
                allocLevels <<= 1;
                Node **alloc = new Node*[allocLevels+1];
                memset(alloc, 0, sizeof(Node*)*(allocLevels+1));
                delete []buckets;
                buckets = alloc;
            }
        }
        inline void free() {
            delete []buckets;
            buckets = NULL;
        }

        template <bool sTree> inline void add(Node* x) {
            int bucket = (sTree ? (x->label) : (-x->label));
            IB_NEXTPTR_EXCESS(x) = buckets[bucket];
            if (buckets[bucket] != NULL) {
                IB_PREVPTR_EXCESS(buckets[bucket]) = x;
            }
            buckets[bucket] = x;
            if (bucket > maxBucket) maxBucket = bucket;
            if (bucket != 0 && bucket < minBucket) minBucket = bucket;
        }
        inline Node* popFront(int bucket) {
            Node *x = buckets[bucket];
            if (x == NULL) return NULL;
            buckets[bucket] = IB_NEXTPTR_EXCESS(x);
            return x;
        }
        template <bool sTree> inline void remove(Node *x) {
            int bucket = (sTree ? (x->label) : (-x->label));
            if (buckets[bucket] == x) {
                buckets[bucket] = IB_NEXTPTR_EXCESS(x);
            } else {
                IB_NEXTPTR_EXCESS(IB_PREVPTR_EXCESS(x)) = IB_NEXTPTR_EXCESS(x);
                if (IB_NEXTPTR_EXCESS(x) != NULL) IB_PREVPTR_EXCESS(IB_NEXTPTR_EXCESS(x)) = IB_PREVPTR_EXCESS(x);
            }
        }
        inline void incMaxBucket(int bucket) {
            if (maxBucket < bucket) maxBucket = bucket;
        }
        inline bool empty() {
            return maxBucket < minBucket;
        }
        inline void reset() {
            maxBucket = 0;
            minBucket = -1 ^ (1<<31);
        }

        Node **buckets;
        Node **ptrs;
        int maxBucket;
        int minBucket;
        Node *nodes;
        int allocLevels;
    };

    Node    *nodes, *nodeEnd;
    Arc *arcs, *arcEnd;
    Node    **ptrs;
    int     numNodes;
    int flow;
    short augTimestamp;
    int topLevelS, topLevelT;
    ActiveList active0, activeS1, activeT1;
    Node **incList;
    int incLen;
    int incIteration;
    Buckets3Pass orphan3PassBuckets;
    BucketsOneSided orphanBuckets;
    ExcessBuckets excessBuckets;
    Buckets3Pass &prNodeBuckets;
    bool verbose;
    
    unsigned int uniqOrphansS, uniqOrphansT;
    template <bool sTree> inline void orphanFree(Node *x) {
        if (IB_EXCESSES && x->excess) {
            x->label = (sTree ? -topLevelT : topLevelS);
            if (sTree) activeT1.add(x);
            else activeS1.add(x);
            x->isParentCurr = 0;
        } else {
            x->label = 0;
        }
    }
    
    struct TmpEdge
    {
        int head;
        int tail;
        int cap;
        int revCap;
    };
    struct TmpArc
    {
        TmpArc *rev;
        int cap;
    };
    char    *memArcs;
    TmpEdge *tmpEdges, *tmpEdgeLast;
    TmpArc *tmpArcs;
    bool isInitializedGraph() {
        return memArcs != NULL;
    }
    void initGraphFast();
    void initGraphCompact();
    void initNodes();
};

inline void MfGraph::addNode(int nodeIndex, int capSource, int capSink)
{
    int f = nodes[nodeIndex].excess;
    if (f > 0) {
        capSource += f;
    } else {
        capSink -= f;
    }
    if (capSource < capSink) {
        flow += capSource;
    } else {
        flow += capSink;
    }
    nodes[nodeIndex].excess = capSource - capSink;
    
}

inline void MfGraph::resetTrees()
{
    resetTrees(1,1);
}

inline void MfGraph::resetTrees(int newTopLevelS, int newTopLevelT)
{
    uniqOrphansS = uniqOrphansT = 0;
    topLevelS = newTopLevelS;
    topLevelT = newTopLevelT;
    for (Node *y=nodes; y != nodeEnd; y++)
    {
        if (y->label < topLevelS && y->label > -topLevelT) continue;
        y->firstSon = NULL;
        if (y->label == topLevelS) activeS1.add(y);
        else if (y->label == -topLevelT) activeT1.add(y);
        else {
            y->parent = NULL;
            if (y->excess == 0) {
                y->label = 0;
            } else if (y->excess > 0) {
                y->label = topLevelS;
                activeS1.add(y);
            } else {
                y->label = -topLevelT;
                activeT1.add(y);
            }
        }
    }
}

inline bool MfGraph::incShouldResetTrees()
{
    return (uniqOrphansS + uniqOrphansT) >= (unsigned int)(2*numNodes);
}

inline void MfGraph::addEdge(int nodeIndexFrom, int nodeIndexTo, int capacity, int reverseCapacity)
{
    tmpEdgeLast->tail = nodeIndexFrom;
    tmpEdgeLast->head = nodeIndexTo;
    tmpEdgeLast->cap = capacity;
    tmpEdgeLast->revCap = reverseCapacity;
    tmpEdgeLast++;
    
    nodes[nodeIndexFrom].label++;
    nodes[nodeIndexTo].label++;
}

inline int MfGraph::isNodeOnSrcSide(int nodeIndex, int freeNodeValue)
{
    if (nodes[nodeIndex].label == 0) {
        return freeNodeValue;
    }
    return (nodes[nodeIndex].label > 0 ? 1 : 0);
}

#endif /* cdn_h */
