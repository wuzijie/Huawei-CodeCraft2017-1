#include "cdn.h"

#define INFO 0                                  // 调试打印模式

#ifdef INFO
#define PRINT   printf
#else
#define PRINT(...)
#endif

int TIME_OUTMIN =  80; // 85s


int totalcost = MAX_NUMBER;       // 当前最大流成本
int totalcost_flow = MAX_NUMBER;  // 最优 方案成本
vector<Edge> g_best_roads;      // 最优 方案路径
vector<ServerNode> g_servers;          // 最优 方案服务器
bitset<MaxNetNode> g_servers_bit,tmp_bit;

int suc_find = 0;

int super_id,join_id,mcmflinkcount,mcmfnodecount; // 超原点超汇点

bool isfindtoplevel;
float netnodeinit_ccc;

TopoInfo topoInfo;      // 图的基本信息
int NetNodeCount;       // 网络节点个数
int ConsumNodeCount;    // 消费节点个数
int LinkCount;          // 链路数
ServerLevel TopServer;  // 最大档次
int TopServerLevel;     // 硬件最大档次
int TopServerFlow;      // 硬件最大档次上限

vector<int> region_root; // 区域父节点索引
float g_slave = 0;


bool servercmp(const ServerNode &a,const ServerNode &b){
    return a.netnode < b.netnode;
}

ServerLevel checkServerLevel(int output){
    for (ServerLevel s : topoInfo.serverLevel) {
        if (output <= s.outPut) {
            return s;
        }
    }
   
    return TopServer;
}

bool TimeOutMin(){
    struct timeb time;
    ftime(&time);
    
    static unsigned long start2 = time.time*1000 + time.millitm;
    unsigned long end = time.time*1000 + time.millitm;
    
    if (end - start2 > TIME_OUTMIN*1000) {
        return true;
    }
    return false;
}

int RunTime(){
    struct timeb time;
    ftime(&time);
    
    static unsigned long start2 = time.time*1000 + time.millitm;
    unsigned long end = time.time*1000 + time.millitm;
    
    return end - start2;
}

void init(){
    
    ConsumNodeCount = topoInfo.ConsumNodeCount;
    NetNodeCount    = topoInfo.NetNodeCount;
    LinkCount       = topoInfo.LinkCount;
    TopServerLevel  = topoInfo.ServerLevelCount-1;
    TopServerFlow = topoInfo.serverLevel[TopServerLevel].outPut;
    TopServer = topoInfo.serverLevel[TopServerLevel];


    super_id = NetNodeCount;
    join_id = NetNodeCount+1;
    mcmflinkcount = LinkCount * 2;
    mcmfnodecount = NetNodeCount + 2;
    
    // 直接连接方案
//    simple();
}

/* 消费节点需求满足情况 */
vector<Edge> tmp_edges; // 本次费用流路由
//int godfailedtime = 0;
//int g_godtime = 0;
bool isScale = false;
bool findRoadsProgrammeWithServerNodebymaxflow(vector<ServerNode> &server) {
//    g_godtime ++;
    tmp_bit.reset();
    totalcost = MAX_NUMBER;
    for (int i = 0; i < server.size(); i++) {
        tmp_bit.set(server[i].netnode);
    }

    // 算上服务器节点
    int num_arcs = mcmflinkcount+(int)server.size()+ConsumNodeCount;
    
    Mcmf mcmf( mcmfnodecount, num_arcs);
    
    int edgecount = (int)topoInfo.g_edges_maxflow.size();
    for (int i = 0; i < edgecount; i++) {
        NetLink &edge = topoInfo.g_edges_maxflow[i];
        mcmf.set_arc(edge.startPoint, edge.endPoint, edge.maxWidth, edge.cost);
        mcmf.set_arc(edge.endPoint, edge.startPoint, edge.maxWidth,  edge.cost);
    }
    
    int a = (int)server.size();
    for (int i = 0; i < a; ++i) { // 超级源。的边
        int level;
        if (isScale==false){
            level = TopServerLevel;
        }else{
            level = min(server[i].serverLevel,TopServerLevel);
        }

        int max;
        if (level == -1)
        {
            max = 0;
        }else{
            max = topoInfo.serverLevel[level].outPut;
        }
        
        mcmf.set_arc(super_id, server[i].netnode, max, 0);
    }
    for (int i = 0; i < ConsumNodeCount; ++i) { // 超汇点。的边
        mcmf.set_arc(topoInfo.consumNode[i].netNode, join_id, topoInfo.consumNode[i].need, 0);
    }
    // 超级汇点 需求
    mcmf.set_supply_demand_of_node(join_id, 0-topoInfo.all_demend);
    // 超级源 供给
    mcmf.set_supply_demand_of_node(super_id,topoInfo.all_demend);
    
    mcmf.exec();
    
    if (mcmf.cost_maxflow >= 0) {
        tmp_edges  = mcmf.roads;
        
        int total= mcmf.cost_maxflow;

        for(vector<ServerNode>::iterator it = server.begin(); it != server.end();it++) {
            if (mcmf.server_supply.find((*it).netnode) != mcmf.server_supply.end()){
                ServerLevel s = checkServerLevel(mcmf.server_supply[(*it).netnode]);
                
                (*it).flow = mcmf.server_supply[(*it).netnode];
                (*it).serverLevel = s.levelid;

                total += (s.cost + topoInfo.net_nodes[(*it).netnode].cost);
                
            }
        }

        totalcost = total;

        if (totalcost_flow > totalcost) {
            g_best_roads = mcmf.roads;
            totalcost_flow = totalcost;
            g_servers = server;
            g_servers_bit = tmp_bit;
            
            suc_find++;

            if (!isfindtoplevel) {
                // 删除多余的服务器
                if (g_servers.size() != mcmf.server_supply.size()){
                    PRINT("DELETE supply:%d server:%d\n",(int)mcmf.server_supply.size(), (int)g_servers.size());
                    for(vector<ServerNode>::iterator it = g_servers.begin(); it != g_servers.end();){
                        if (mcmf.server_supply.find((*it).netnode) == mcmf.server_supply.end()){
                            it = g_servers.erase(it);
                        }
                        else{
                            it++;
                        }
                    }
                }
            }
            
        }
        
        return true;
    } else {
//        godfailedtime++;
        return false;
    }
}

//------------------------------------------------------------
int pain = 0;
int randomi(int a, int b)
{
    int c=(rand()+pain)%(b-a+1)+a;
    return c;
}
double randomf(double a, double b)
{
    
    double c = (double)(rand()%((int)b-(int)a)) + a + (double)(rand()/(RAND_MAX + 1.0));
    return c;
}

class Solution{
public:
    friend bool operator < (Solution n1, Solution n2){
        return n1.cost > n2.cost;
    }
    Solution(vector<ServerNode> s,int c) {
        server = s;
        cost = c;
    }
    Solution(vector<ServerNode> s,vector<Edge> r,int c) {
        server = s;
        roads = r;
        cost = c;
    }
    Solution(vector<ServerNode> s,vector<Edge> r,int c,int l) {
        server = s;
        roads = r;
        cost = c;
        toplevel = l;
    }
    Solution(vector<ServerNode> s,vector<Edge> r,int c,bitset<MaxNetNode> b) {
        server = s;
        roads = r;
        cost = c;
        bit = b;
    }
    Solution() {
        cost = MAX_NUMBER;
    }
    
    bool runMaxflow() {
        bool suc = findRoadsProgrammeWithServerNodebymaxflow(server);
        if (suc == false || totalcost > totalcost_flow * g_slave) {
            return false;
        } else {
            roads   = tmp_edges;
            bit     = tmp_bit;
            cost    = totalcost;
            return true;
        }
    }
    vector<ServerNode> server;
    vector<Edge> roads;
    int cost;
    bitset<MaxNetNode> bit;
    int toplevel;
};


RouteDetail compute_unit_cost(vector<Edge> &route_edge, bitset<MaxNetNode> &server_bit, int type);

bool findMaxFlow(vector<ServerNode> &server);

void xjbdel() {
    // 删除
    isScale = true;
    int nnum;//网络节点个数
    vector<ServerNode> server = compute_unit_cost(g_best_roads, g_servers_bit, 1).server_set;
    int a = (int)server.size();
    for (int i = 0; i < a; ++i)
    {
        if (topoInfo.key_bit.test(server[i].netnode)){
           nnum = i;
           break;
        }
    } 
 
    int cnum = (int)server.size() - nnum; // 消费节点个数
    
    int stop_num = (int)nnum * 0.71; // 删网络节点个数，
    int next_start = nnum - stop_num;// 删消费节点起点，
   
    while(stop_num && !TimeOutMin()) {
        ServerNode temp = server[0];
        server.erase(server.begin());

        if (findMaxFlow(server)){
            findRoadsProgrammeWithServerNodebymaxflow(server);
            if (totalcost > totalcost_flow) { // 没有更新最优
                server.push_back(temp);
            }else{
                PRINT("better del %d totalcost: %d\n", temp.netnode, totalcost_flow);
            }
        }
        else{ // 删除失败
            server.push_back(temp);
        }

        stop_num--;
    }

    stop_num = (int) cnum * 0.37;
    while(stop_num && !TimeOutMin()) {
        ServerNode temp = server[next_start];
        server.erase(server.begin()+next_start);

        if (findMaxFlow(server)){
            findRoadsProgrammeWithServerNodebymaxflow(server);
            if (totalcost > totalcost_flow) { // 没有更新最优
                server.push_back(temp);
            }else{
                PRINT("better del %d totalcost: %d\n", temp.netnode, totalcost_flow);
            }
        }
        else{ // 删除失败
            server.push_back(temp);
        }
        
        stop_num--;
    }
}


bool swapNode_d_f_c(const ServerNode &a, const ServerNode &b){ // 给网络节点排序，，先度 后输出能力 再成本
    NetNode a_ = topoInfo.net_nodes[a.netnode];
    NetNode b_ = topoInfo.net_nodes[b.netnode];
    
    if (a_.flow_ ==  b_.flow_) {
        if (a_.degree_ == b_.degree_) {
            return a_.cost < b_.cost;
        }
        return a_.degree_ > b_.degree_;
    }
    return a_.flow_ > b_.flow_;
}


// 排序
bool comp_xjbs(const ServerNode &a, const ServerNode &b)
{
    if (topoInfo.net_nodes[a.netnode].flow_ == topoInfo.net_nodes[b.netnode].flow_){
        return topoInfo.net_nodes[a.netnode].degree_ > topoInfo.net_nodes[b.netnode].degree_;
    }
    return topoInfo.net_nodes[a.netnode].flow_ > topoInfo.net_nodes[b.netnode].flow_;
}


void xjbs() {
    int sideNum = 50;
    if (NetNodeCount > 500)
        sideNum = 25;
    if (NetNodeCount > 1000)
        sideNum = 60;
    
    vector<ServerNode> server = g_servers;
    sort(server.begin(),server.end(),comp_xjbs);
    
    
    Solution s(server,g_best_roads,totalcost_flow,g_servers_bit);
    priority_queue<Solution> q;
    q.push(s);
    for (int j = 0; j < s.server.size() && !TimeOutMin(); j++) {
        Solution s = q.top();
        priority_queue<Solution> back;
        int index = j;
        // iter sidenode
        NetLink *link  = topoInfo.graph[s.server[index].netnode];
        vector<ServerNode> swapNode;
        while (link!=NULL) {
            int node = link -> endPoint;
            link = link -> next;
            if (s.bit.test(node)) {
                continue;
            }
            swapNode.push_back(ServerNode(node, topoInfo.net_nodes[node].maxServerLevel));
            
        }
        
        sort(swapNode.begin(), swapNode.end(), swapNode_d_f_c);
        
        for (int k = 0; k <min(int(swapNode.size()),sideNum); k++) {
            int index2 = k;
            Solution temp = s;
            swap(temp.server[index].netnode,swapNode[index2].netnode);
            if (temp.runMaxflow()) {
                PRINT("betterswap %d to %d = %d \n",temp.server[index].netnode,swapNode[index2].netnode,totalcost);
                back.push(temp);
            } else {
                swap(temp.server[index].netnode,swapNode[index2].netnode);
            }
        }
        if (back.empty()) back.push(s);
        q = back;
    }
}


void xjbs(int step) {
    int sideNum = 50;
    if (NetNodeCount > 500)
        sideNum = 60;
    if (NetNodeCount > 1000)
        sideNum = 80;
    
    PRINT("\ncurrentnodesize: %d cost %d \n",(int)g_servers.size(),totalcost_flow);
    
    Solution s(g_servers,g_best_roads,totalcost_flow,g_servers_bit);
    
    priority_queue<Solution> q;
        
    RouteDetail detail = compute_unit_cost(s.roads, s.bit, 1);
    vector<ServerNode> needSwapNode = detail.server_set;
    

//    int maxcount = (int)(needSwapNode.size()*(0.91 - 0.05*step));
    int maxcount = (int)needSwapNode.size();
    s.server = needSwapNode;
    q.push(s);
    for (int i = 0; i < maxcount && !TimeOutMin(); i++) {
        Solution s = q.top();
        priority_queue<Solution> back;
        
        if (!s.bit.test(needSwapNode[i].netnode)) {
            continue;
        }
        int index = i;
        // iter sidenode
        NetLink *link  = topoInfo.graph[s.server[index].netnode];
        vector<ServerNode> swapNode;
        while (link!=NULL) {
            int node = link -> endPoint;
            link = link -> next;
            if (s.bit.test(node)) {
                continue;
            }
            swapNode.push_back(ServerNode(node, topoInfo.net_nodes[node].maxServerLevel));
        }
        
        sort(swapNode.begin(), swapNode.end(), swapNode_d_f_c);
        
        for (int k = 0; k <min(int(swapNode.size()),sideNum); k++) {
            int index2 = k;
            Solution temp = s;
            swap(temp.server[index].netnode,swapNode[index2].netnode);
            if (temp.runMaxflow()) {
                PRINT("betterswap %d to %d = %d \n",temp.server[index].netnode,swapNode[index2].netnode,totalcost);
                back.push(temp);
            }
        }
        
        if (back.empty())
            back.push(s);
        
        q = back;
    }
    
}

void xjbs_route(int step) {
    int sideNum = 50;
    if (NetNodeCount > 500)
        sideNum = 25;
    if (NetNodeCount > 1000)
        sideNum = 60;
    
    PRINT("\ncurrentnodesize: %d cost %d \n",(int)g_servers.size(),totalcost_flow);
    
    Solution s(g_servers,g_best_roads,totalcost_flow,g_servers_bit);
    
    priority_queue<Solution> q;
        
    RouteDetail detail = compute_unit_cost(s.roads, s.bit, 4);
    vector<ServerNode> needSwapNode = detail.server_set;
    unordered_map<int, unordered_set<int>> route_set = detail.route_set; 
    

//    int maxcount = (int)(needSwapNode.size()*(0.91 - 0.05*step));
    int maxcount = (int)needSwapNode.size();
    s.server = needSwapNode;
    q.push(s);

    for (int i = 0; i < maxcount && !TimeOutMin(); i++) {
        Solution s = q.top();
        priority_queue<Solution> back;
        
        if (!s.bit.test(needSwapNode[i].netnode)) {
            continue;
        }

        int index = i;
        // iter sidenode
       
        vector<ServerNode> swapNode;
        for (int node : route_set[needSwapNode[i].netnode]) {
            if (s.bit.test(node)) {
                continue;
            }
            swapNode.push_back(ServerNode(node, topoInfo.net_nodes[node].maxServerLevel));
        }
        
        sort(swapNode.begin(), swapNode.end(), swapNode_d_f_c);
        
        for (int k = 0; k <min(int(swapNode.size()),sideNum); k++) {
            int index2 = k;
            Solution temp = s;
            swap(temp.server[index].netnode,swapNode[index2].netnode);
            if (temp.runMaxflow()) {
                PRINT("route betterswap %d to %d = %d \n",temp.server[index].netnode,swapNode[index2].netnode,totalcost);
                back.push(temp);
            }
        }
        
        if (back.empty())
            back.push(s);
        
        q = back;
    }
    
}

void xjbdownlevel(){
    // 降档次
    // 降最单位成本高的档次， （单位成本高：1.消费节点 2.没充分利用档次）
    PRINT("-----------xjbdownlevel-----------\n");
    
    isScale =true;
    
    vector<ServerNode> server = compute_unit_cost(g_best_roads, g_servers_bit, 1).server_set;
    
    int net_stop = 0;// consum_start,consum_stop;
    
    int a = (int)server.size();
    for (int i = 0; i < a && !TimeOutMin(); i++) {
        if (topoInfo.key_bit.test(server[i].netnode)) {
            net_stop = i+1;
            break;
        }
    }
    
    net_stop *= 0.87;

    for (int i = 0 ; i < net_stop && !TimeOutMin();i++) {
        
        int level = server[i].serverLevel--;
        if (findMaxFlow(server)) {
            findRoadsProgrammeWithServerNodebymaxflow(server);
            if (totalcost > totalcost_flow) {
                server[i].serverLevel++;//
            } else {
                
                PRINT("downlevel server%d  %d to %d = %d \n",server[i].netnode,level,level-1,totalcost);
            }
        }else{
            server[i].serverLevel++;
        }
        
    }
    
    PRINT("----------------------------\n");
}

void xjbuplevel(){
    // 升档次
    PRINT("-----------xjbuplevel-----------\n");
    
    isScale =true;
    vector<ServerNode> server = g_servers;
    
    bitset<MaxNetNode> up_bit;
    for (int v : compute_unit_cost(g_best_roads, g_servers_bit, 3).up_set){
        up_bit.set(v);
    }
          
    PRINT("totalcost_flow: %d\n", totalcost_flow);
    int a = (int)server.size();
    for (int i = 0;i < a && !TimeOutMin();i++){
        if (up_bit.test(server[i].netnode)){
            server[i].serverLevel++;
            bool isSuc = findRoadsProgrammeWithServerNodebymaxflow(server);
            if (totalcost > totalcost_flow || isSuc == false) {// 删除失败
                server[i].serverLevel--;
            }else{
                PRINT("better up %d totalcost_flow: %d\n", server[i].netnode, totalcost_flow);
            }
        }
    }
    
    PRINT("----------------------------\n");
}

void writeresult(char * filename){
    // 输出
    string result;
    
    PRINT("1 server:");
    for (int i = 0; i < g_servers.size(); i++) {
       PRINT("%d, ",g_servers[i].netnode);
    }
    PRINT("\n");
    PRINT("2 cost:%d\n",totalcost_flow);
    PRINT("3 level:");
    for (int i = 0; i < g_servers.size(); i++) {
       PRINT("%d, ",g_servers[i].serverLevel);
    }
    PRINT("\n");
    
    
    writeMaxFlowPrograme(result, g_servers, g_best_roads, NetNodeCount, topoInfo.g_net2consume);

    char * topo_file = (char *)result.c_str();
    write_result(topo_file, filename);
}


/************************************* 最大流 **********************/
vector<int> fail_consumer;
void analysisFailMaxFlow(vector<ServerNode> &server, int real_flow);
bool analysis_fail_flow = false;

/**
 * [findMaxFlow 计算服务器方案最大流]
 * @return [是否可行]
 */
bool findMaxFlow(vector<ServerNode> &server){
    //算上超级源点
//    int num_nodes = NetNodeCount + 2;// ljx
    // 算上服务器节点
    int num_arcs = mcmflinkcount + (int)server.size() + ConsumNodeCount; // ljx
    
    MfGraph mf;
    mf.initSize(mcmfnodecount, num_arcs);
    for (int i = 0; i < NetNodeCount; i++){
        mf.addNode(i, 0, 0);
    }
    mf.addNode(super_id,topoInfo.all_demend, 0);
    mf.addNode(join_id, 0, topoInfo.all_demend);
    
    int edgecount = (int)topoInfo.g_edges_maxflow.size();
    for (int i = 0; i < edgecount; i++) {
        NetLink &edge = topoInfo.g_edges_maxflow[i];
        mf.addEdge(edge.startPoint, edge.endPoint, edge.maxWidth, 0);
        mf.addEdge(edge.endPoint, edge.startPoint, edge.maxWidth, 0);
    }
    TopServerFlow = topoInfo.serverLevel[TopServerLevel].outPut;
    TopServer = topoInfo.serverLevel[TopServerLevel];
    for (ServerNode s : server) { // 超级源。的边
        if (isScale)
            mf.addEdge(super_id, s.netnode, topoInfo.serverLevel[s.serverLevel].outPut, 0);
        else
            mf.addEdge(super_id, s.netnode, TopServerFlow, 0);
    }
    for (int i = 0; i < ConsumNodeCount; ++i) { // 超汇点。的边
        mf.addEdge(topoInfo.consumNode[i].netNode, join_id, topoInfo.consumNode[i].need, 0);
    }
    mf.initGraph();
    mf.computeMaxFlow();
    
    int max_flow = mf.getFlow();
    if (max_flow != topoInfo.all_demend){
        if (analysis_fail_flow)  analysisFailMaxFlow(server, max_flow);
        return false;
    } 
    else return true;
}

void analysisFailMaxFlow(vector<ServerNode> &server, int real_flow){
    //算上超级源点
//    int num_nodes = NetNodeCount + 2;// ljx
    // 算上服务器节点
    int num_arcs = mcmflinkcount + (int)server.size() + ConsumNodeCount; // ljx
    
    MfGraph mf;
    mf.initSize(mcmfnodecount, num_arcs);

    for (int i = 0; i < NetNodeCount; i++){
        mf.addNode(i, 0, 0);
    }
    mf.addNode(super_id, real_flow, 0);
    mf.addNode(join_id, 0, real_flow);

    int edgecount = (int)topoInfo.g_edges_maxflow.size();
    for (int i = 0; i < edgecount; i++) {
        NetLink &edge = topoInfo.g_edges_maxflow[i];
        mf.addEdge(edge.startPoint, edge.endPoint, edge.maxWidth, 0);
        mf.addEdge(edge.endPoint, edge.startPoint, edge.maxWidth, 0);
    }
    for (ServerNode s : server) { // 超级源。的边
        if (isScale)
            mf.addEdge(super_id, s.netnode, topoInfo.serverLevel[s.serverLevel].outPut, 0);
        else
            mf.addEdge(super_id, s.netnode, TopServerFlow, 0);;
    }
    for (int i = 0; i < ConsumNodeCount; ++i) { // 超汇点。的边
        mf.addEdge(topoInfo.consumNode[i].netNode, join_id, topoInfo.consumNode[i].need, 0);
    }

    mf.initGraph();
    mf.computeMaxFlow();
    mf.analysis_fail();
    fail_consumer = mf.fail_consumer;
}

/************************************* 最大流 **********************/


//----type 1234 排序
int cmp_unit_net(const pair<int, float> &x, const pair<int, float> &y){
    return x.second < y.second;
}
int cmp_unit_server(const pair<int, float> &x, const pair<int, float> &y){
    if (!topoInfo.key_bit.test(x.first) && topoInfo.key_bit.test(y.first)) {
        return true;
    }
    if (topoInfo.key_bit.test(x.first) && !topoInfo.key_bit.test(y.first)) {
        return false;
    }
    
    return x.second > y.second;
}
int cmp_unit_pass(const pair<int, int> &x, const pair<int, int> &y){
    return x.second > y.second;
}


/**
 * 分析费用流
 * route_edge 费用流的路由边
 */
RouteDetail compute_unit_cost(vector<Edge> &route_edge, bitset<MaxNetNode> &server_bit, int type) {
    // 存储路径
    unordered_map<int, vector<Edge> > road_map;
    unordered_map<int, Edge> other_road_map;
    unordered_map<int, int> net_cost, server_cost;
    unordered_map<int, int> net_width, server_width;
    unordered_map<int, int> pass_server;    // 有流量流经的 服务器集
    unordered_map<int, unordered_set<int>> route_set;    // 服务器的路由节点集
    unordered_map<int, unordered_set<int>> service_set;    //  服务器服务的消费节点集 

    RouteDetail route_detail;
    route_detail.type = type;
    
    for (Edge e : route_edge) {
        if (road_map.find(e.start) == road_map.end())  // 寻找路径用
            road_map[e.start] = vector<Edge>();
        
        road_map[e.start].push_back(e);
        other_road_map[e.start*10000+e.end] = e; // 分配路径用
    }
    queue<Flow> q;
    Flow f;
    f.cur_node = super_id; // NetNodeCount 为超级源点编号
    f.roads.push_back(super_id);

    q.push(f); // 超级源
    vector<vector<int>> roads;  // 寻找出来的路径  最后一个是超级汇点  下标与 answers。对应

    while (!q.empty()) {
        Flow f = q.front();
        q.pop();
        
        if (f.cur_node == join_id) { // 到超汇点
            roads.push_back(f.roads);   //记录 路径 用于分配带宽
            continue;
        }
        
        vector<Edge> v = road_map[f.cur_node];
        
        for (Edge e : v) {
            Flow tmp = f;
            tmp.roads.push_back(e.end);
            tmp.cur_node = e.end;

            q.push(tmp);
        }
    }
    
    for (vector<int> r : roads) {  /// 在这里确定 路径的 宽带。。 路径上的链路有些为0的时候就删除这条路
        int minWidth = MAX_NUMBER;
        
        for (int j = 0; j < r.size()-1; j++) {
            int start   = r[j];
            int end     = r[j+1];
           
            if (other_road_map[start*10000+end].flow < minWidth) { // 找最小
                minWidth = other_road_map[start*10000+end].flow;
            }
        }

        if (minWidth != 0) { // 可以分配 带宽
            unordered_map<int, int> tmp_cost;   
            int cost = 0;
            int sid = r[1], cid = r[r.size()-2];

            for (int j = 0; j < r.size()-1; j++) { // 减去。使用过的
                int start   = r[j];
                int end     = r[j+1];
                other_road_map[start*10000+end].flow -= minWidth;
                
                if (j > 1 && j < r.size()-2) {
                    if (!server_bit.test(start)){
                        tmp_cost[start] = cost;

                        if(type == 4) {
                            route_set[sid].insert(start);
                        }
                    }
                    else {
                        if (type == 3){
                            if (pass_server.find(start) == pass_server.end()){
                                pass_server[start] = minWidth;
                            }
                            else{
                                pass_server[start] += minWidth;
                            }
                        }
                    }
                }

                cost += other_road_map[start*10000+end].rent;
            }

            for (auto it : tmp_cost){
                if (net_cost.find(it.first) == net_cost.end()){
                    net_cost[it.first] = (cost - it.second) * minWidth;
                    net_width[it.first] = minWidth;
                }
                else{
                    net_cost[it.first] += (cost - it.second) * minWidth;
                    net_width[it.first] += minWidth;
                }
            }

            
            if (server_cost.find(sid) == server_cost.end()){
                server_cost[sid] = cost * minWidth;
                server_width[sid] = minWidth;
            }
            else{
                server_cost[sid] += cost * minWidth;
                server_width[sid] += minWidth;
            }

            service_set[sid].insert(cid);
        }
    }

    vector<pair<int, float>> pass_set, server_set;
    for (auto it : net_cost){
        float width_unit;
        if (topoInfo.key_bit.test(it.first)){
            int cid = topoInfo.g_net2consume[it.first];
            int width_consumer = net_width[it.first]+topoInfo.consumNode[cid].need;

            if (width_consumer > TopServerFlow){
                int all_cost = it.second + topoInfo.net_nodes[it.first].cost + TopServer.cost;
                width_unit = all_cost / (float) (TopServerFlow);
            }
            else {
                int all_cost = it.second + topoInfo.net_nodes[it.first].cost + checkServerLevel(width_consumer).cost;
                width_unit = all_cost / (float) (width_consumer);
            }
        }
        else{
            if (net_width[it.first] > TopServerFlow) {
                int all_cost = it.second + topoInfo.net_nodes[it.first].cost + TopServer.cost;
                width_unit = all_cost / (float) (TopServerFlow);
            }
            else {
                int all_cost = it.second + topoInfo.net_nodes[it.first].cost + checkServerLevel(net_width[it.first]).cost;
                width_unit = all_cost / (float) (net_width[it.first]);
            }
        }

        pass_set.push_back(pair<int, float>(it.first, width_unit));
    }

    for (auto it : server_cost){
        int all_cost = it.second + topoInfo.net_nodes[it.first].cost + checkServerLevel(server_width[it.first]).cost;
        float width_unit = all_cost / (float) server_width[it.first];

        server_set.push_back(pair<int, float>(it.first, width_unit));
    }

    if (type == 1)
        sort(server_set.begin(), server_set.end(), cmp_unit_server);
    else if (type == 2)
        sort(pass_set.begin(), pass_set.end(), cmp_unit_net);
    else if (type == 3){
        vector<pair<int, int>> pass_pair(pass_server.begin(), pass_server.end());
       
        for(auto it : pass_pair){
            if (checkServerLevel(server_width[it.first]).levelid < TopServerLevel) {
                route_detail.up_set.push_back(it.first);
            }
        }
    }
    else{
        route_detail.route_set = route_set;
        sort(server_set.begin(), server_set.end(), cmp_unit_server);
    }

    
    for (auto it : server_set){
        topoInfo.net_nodes[it.first].unit_cost = it.second;
        if (type == 1 || type == 4){
            route_detail.server_set.push_back(ServerNode(it.first, checkServerLevel(server_width[it.first]).levelid, server_width[it.first]));
        }

        // cout<<"# "<<it.first<<" unit: "<<it.second<<" cost: "<<topoInfo.net_nodes[it.first].cost<<" level: "<<checkServerLevel(server_width[it.first]).levelid
        //         <<" route: "<<server_cost[it.first]<<" width: "<<server_width[it.first]<<" degree: "<<topoInfo.net_nodes[it.first].degree_<<"  flow: "<<topoInfo.net_nodes[it.first].flow_<<endl;
    }

    for (auto it : pass_set){
        topoInfo.net_nodes[it.first].unit_cost = it.second;
        if (type == 2){
            if (!topoInfo.key_bit.test(it.first)) {
                if (net_width[it.first] > TopServerFlow)
                    route_detail.pass_set.push_back(ServerNode(it.first, TopServerLevel, TopServerFlow));
                else
                    route_detail.pass_set.push_back(ServerNode(it.first, checkServerLevel(net_width[it.first]).levelid, net_width[it.first]));
            }
            else{
                int cid = topoInfo.g_net2consume[it.first];
                int width_consumer = net_width[it.first]+topoInfo.consumNode[cid].need;

                if (width_consumer > TopServerFlow){
                    route_detail.pass_set.push_back(ServerNode(it.first, TopServerLevel, TopServerFlow));
                }
                else {
                    route_detail.pass_set.push_back(ServerNode(it.first, checkServerLevel(width_consumer).levelid, width_consumer));
                }
            }
        }

        // cout<<"* "<<topoInfo.key_bit.test(it.first)<<" "<<it.first<<" unit: "<<it.second<<" cost: "<<topoInfo.net_nodes[it.first].cost<<" level: "<<checkServerLevel(net_width[it.first]).levelid
        //     <<" route: "<<net_cost[it.first]<<" width: "<<net_width[it.first]<<" degree: "<<topoInfo.net_nodes[it.first].degree_<<"  flow: "<<topoInfo.net_nodes[it.first].flow_<<endl;
    }

    return route_detail;
}


/************************** 二进制初始解 *****************************/
bool comp_binary_degree(const ServerNode &a, const ServerNode &b)
{
    if (topoInfo.net_nodes[a.netnode].degree_ == topoInfo.net_nodes[b.netnode].degree_){
        if (topoInfo.net_nodes[a.netnode].cost == topoInfo.net_nodes[b.netnode].cost){
            return topoInfo.net_nodes[a.netnode].flow_ < topoInfo.net_nodes[b.netnode].flow_;
        }
        return topoInfo.net_nodes[a.netnode].cost > topoInfo.net_nodes[b.netnode].cost;
    }
    return topoInfo.net_nodes[a.netnode].degree_ < topoInfo.net_nodes[b.netnode].degree_;
}


bool comp_binary_flow(const ServerNode &a, const ServerNode &b)
{
    if (topoInfo.net_nodes[a.netnode].flow_ == topoInfo.net_nodes[b.netnode].flow_){
         if (topoInfo.net_nodes[a.netnode].cost == topoInfo.net_nodes[b.netnode].cost){
            return topoInfo.net_nodes[a.netnode].degree_ < topoInfo.net_nodes[b.netnode].degree_;
        }
        return topoInfo.net_nodes[a.netnode].cost > topoInfo.net_nodes[b.netnode].cost;
    }
    return topoInfo.net_nodes[a.netnode].flow_ < topoInfo.net_nodes[b.netnode].flow_;
}

bool comp_binary_unit(const ServerNode &a, const ServerNode &b)
{
    if (topoInfo.net_nodes[a.netnode].unit_cost == topoInfo.net_nodes[b.netnode].unit_cost){
        if (topoInfo.net_nodes[a.netnode].degree_ == topoInfo.net_nodes[b.netnode].degree_){
            return topoInfo.net_nodes[a.netnode].flow_ < topoInfo.net_nodes[b.netnode].flow_;
        }
        return topoInfo.net_nodes[a.netnode].degree_ < topoInfo.net_nodes[b.netnode].degree_;
    }
    return topoInfo.net_nodes[a.netnode].unit_cost > topoInfo.net_nodes[b.netnode].unit_cost;
}

void binaryinit(int type, vector<ServerNode> init_server) {
    totalcost_flow = MAX_NUMBER;
    vector<ServerNode> server,tmp,tmp2,back;
    server.assign(init_server.begin(), init_server.end());

    if (type == 1)
        sort(server.begin(),server.end(),comp_binary_degree);
    else if (type == 2)
        sort(server.begin(),server.end(),comp_binary_flow);
    else
        sort(server.begin(),server.end(),comp_binary_unit);

    int num = 1;
    
    PRINT("type:%d del start with size %d \n", type, (int)server.size());
    while (true) {
        tmp.clear();
        for (int i = num; i < server.size(); i++) {
            tmp.push_back(server[i]);
        }
        tmp2 = tmp;

        for (ServerNode s : back)
            tmp2.push_back(s);

        //PRINT("size: %d num:%d back:%d tmp:%d tmp2:%d is_feasible:%d\n", server.size(), num, back.size(), tmp.size(), tmp2.size(), is_feasible);
        if (findMaxFlow(tmp2)) {
            //PRINT("del suc with size : %d \n",tmp2.size());
            server = tmp;
            num *= 2;
        } else {
            if (num != 1) {
                num = 1;
            } else {
                back.push_back(server.front());
            }
        }
        
        if (server.empty()) break;
    }
   
    if (findRoadsProgrammeWithServerNodebymaxflow(tmp2)) {
        PRINT("after del size:%d cost : %d cost: %d\n", (int)tmp2.size(), totalcost, totalcost_flow);
    }
}

void netnodeinit(Solution solu, int step){
    totalcost = MAX_NUMBER;
    totalcost_flow = solu.cost;
    g_servers = solu.server;
    g_servers_bit = solu.bit;
    float reduce;

    if (step == 0){
        reduce = 1.0;
    }
    else{
//        reduce = pow(0.91,step);
         reduce = 0.91 - 0.07*step;
    }
    
    int accpet_offset = totalcost_flow * netnodeinit_ccc * reduce;
    
    vector<ServerNode> server;
    PRINT("step: %d reduce:%f init size:%d cost: %d\n", step, reduce, (int)solu.server.size(), solu.cost);
    server = solu.server;
   
    RouteDetail route_detail = compute_unit_cost(solu.roads, g_servers_bit, 2);
    vector<ServerNode> choosenode = route_detail.pass_set;

    if (step == 0) {
        /* 去掉度为2的网络节点 */
        for(vector<ServerNode>::iterator it = choosenode.begin(); it != choosenode.end();) {
            if (!topoInfo.key_bit.test((*it).netnode) && topoInfo.net_nodes[(*it).netnode].degree_ == 2){
                it = choosenode.erase(it);
            }
            else
                it++;
        }
    }
   
    int iter_num = (int) choosenode.size() * 0.67 * reduce;

    //  添加
    for (int i = 0; i < iter_num && !TimeOutMin();i++) {
        if (step != 0) {
            if (topoInfo.key_bit.test(choosenode[i].netnode)){
                if (iter_num != choosenode.size())
                    iter_num++;
                continue;
            }
        }
        
        choosenode[i].serverLevel = TopServerLevel;
        server.push_back(choosenode[i]);
        bool isSuc = findRoadsProgrammeWithServerNodebymaxflow(server);
        if (isSuc == false || totalcost > totalcost_flow + accpet_offset) {// add失败
            server.pop_back();
        }
        else{
            if(totalcost != totalcost_flow){
                g_best_roads = tmp_edges;
                totalcost_flow = totalcost;
                g_servers = server;
                g_servers_bit = tmp_bit;
            }
        }
    }
    
    PRINT("serversize:%d\t g_size:%d\t totalcost:%d\n", (int)server.size(), (int)g_servers.size(), totalcost_flow);
}



bool fireAccept(int deta /*bestcost - totalcost*/,float &T) {
    // 1. if xxx return true;
    // 2. else
    
    float tips = 0.0;// ***
    
    if (deta > 0){
        return true;
    }
    else if(((deta<0)&&(exp(deta/T)-tips>randomf(0,1)))){
        return true;
    }
    return false;
}
void fireScaleDown(vector<ServerNode> &initServers,int &min,float &T) {
    // rand scaleup one from initServers
    // 1. rand
    // 2. maxflow -> false or true
    // 3. maxflow == true ->  maxflowmincost false or true
    
    int index = randomi(0,(int)initServers.size()-1);
    if (initServers[index].serverLevel == 0) {
        return;
    }
    //PRINT("%d scale down %d\n",index,min);
    initServers[index].serverLevel--; // down
    
    if (findMaxFlow(initServers)) {
        min = totalcost_flow;
        int suc = findRoadsProgrammeWithServerNodebymaxflow(initServers);
        if (fireAccept(min-totalcost,T)  && suc == true) {   // suc
            PRINT("fireScaleDown cost : %d\n",totalcost);
            min = totalcost;
        } else {                                        // failed
            initServers[index].serverLevel++;
        }
    } else {
        initServers[index].serverLevel++;
    }
    
}
void fireScaleUp(vector<ServerNode> &initServers,int &min,float &T) {
    // rand scaleup one from initServers
    // 1. rand
    // 2. maxflowmincost -> false or true
    
    // 1. rand
    int index = randomi(0,(int)initServers.size()-1);
    if (initServers[index].serverLevel == TopServerLevel) {
        return;
    }
    //PRINT("scale up\n");
    initServers[index].serverLevel++; // up
    min = totalcost_flow;
    int suc = findRoadsProgrammeWithServerNodebymaxflow(initServers);
    if (fireAccept(min-totalcost,T)  && suc == true) {   // suc
        PRINT("fireScaleUp cost : %d\n",totalcost);
        min = totalcost;
    } else {                                        // failed
        initServers[index].serverLevel--;
    }
}
void fireCore(vector<ServerNode> tempInit,vector<ServerNode> tempChoose) {
    // init answers;
    vector<ServerNode> initServers = tempInit;
    vector<ServerNode> chooseServers = tempChoose;
    
    // init params;
    int min = totalcost_flow;
    float T = 100, f = 0.985, Limit = 0.0000001;
    if (NetNodeCount > 1000) {
        T = 1000, f = 0.98, Limit = 0.00000001;
    }
    int step = 0;
    // iterator
    PRINT("Start fire\n");
    while (T > Limit&&!TimeOut()) {
        //        PRINT("initservers size %d\n", (int)initServers.size());
        step++;
        int index = randomi(0,7);
        if (index <= 3) {
            fireScaleUp(initServers,min,T);
        } else {
            fireScaleDown(initServers,min,T);
        }
        T *= f;
    }
    PRINT("End fire\n");
}

priority_queue<Solution> q;
int iter_size = 0;
void middleCase() {
    isScale = false;
    vector<ServerNode> init_server;
    for (int i = 0; i < NetNodeCount; i++) { // 直接链接的服务器
        if (!topoInfo.key_bit.test(i) && topoInfo.net_nodes[i].degree_ < 3) continue; // wh
        init_server.push_back(ServerNode(i, TopServerLevel));
    }

    unordered_map<int, Solution> init_solu;
    for (int i = 1; i < 4; ++i)
    {
        binaryinit(i, init_server);
        Solution s(g_servers, g_best_roads, totalcost_flow, g_servers_bit);
        init_solu[i] = s;
    }

    for (int i = 1; i < 4; i++) {
        netnodeinit(init_solu[i], 0);
        q.push(Solution(g_servers, g_best_roads, totalcost_flow, g_servers_bit));
    }
    //return;
    isScale = true;

    priority_queue<Solution> p = q;
    q = priority_queue<Solution>();
    int prevCost = MAX_NUMBER;
    while (!TimeOutMin() && p.empty()==false) {
        iter_size++;
        Solution st = p.top();
        p.pop();
        if (st.cost == prevCost) continue;
        prevCost = st.cost;
        g_servers = st.server;
        g_best_roads = st.roads;
        totalcost_flow = st.cost;
        g_servers_bit = st.bit;
        
        for (int i = 0; i < 8 && !TimeOutMin(); i++) {
            PRINT("-----------------step %d-----------------\n", i+1);
            if (i == 0) {
                isScale = false;
                xjbs(0);
                xjbs_route(0);
                isScale = true;
                xjbdel();
                fireCore(g_servers, vector<ServerNode>());
            }
            //if (i == 0) {
                isScale = false;
                bool ok = true;
                while (ok) {
                    int tmin = totalcost_flow;
                    int tlevel = TopServerLevel;
                    TopServerLevel--;
                    bool suc = findRoadsProgrammeWithServerNodebymaxflow(g_servers);
                    if (suc == false || totalcost > tmin) {
                        TopServerLevel = tlevel;
                        ok = false;
                        break;
                    }
                }
            //}
            if (TimeOutMin()) break;
            if (i % 2 == 0)
                xjbs(i);
            else
                xjbs_route(i);
            
            xjbdel();
            if (TimeOutMin()) break;
        }
        //------------
        
        q.push(Solution(g_servers, g_best_roads, totalcost_flow,g_servers_bit));
        PRINT("*******totalcost:%d********\n",totalcost_flow);
    }
}

void change(){
    int all_cost, flow_limit;
    for (int i = 0; i < NetNodeCount; i++) {
        all_cost = 0;
        flow_limit = 0;
        
        int output = topoInfo.net_nodes[i].flow_;
        if(output >= TopServerFlow){
            topoInfo.net_nodes[i].maxServerLevel = TopServerLevel;
            all_cost += TopServer.cost;
            flow_limit = TopServerFlow;
        } else {
            ServerLevel s = checkServerLevel(output);
            topoInfo.net_nodes[i].maxServerLevel = s.levelid;
            all_cost += s.cost;
            flow_limit = s.outPut;
        }
        
        all_cost += topoInfo.net_nodes[i].cost;
        
        if (topoInfo.key_bit.test(i)) {
            int cid = topoInfo.g_net2consume[i];
            if (topoInfo.net_nodes[i].flow_ <= flow_limit){
                all_cost += (topoInfo.net_nodes[i].flow_ - topoInfo.consumNode[cid].need) * topoInfo.width_unit;
                topoInfo.net_nodes[i].unit_cost = all_cost / (float) topoInfo.net_nodes[i].flow_;
            }
            else{
                if(topoInfo.consumNode[cid].need < flow_limit){
                    all_cost += (flow_limit - topoInfo.consumNode[cid].need) * topoInfo.width_unit;
                    topoInfo.net_nodes[i].unit_cost = all_cost / (float) flow_limit;
                }
                else{
                    topoInfo.net_nodes[i].unit_cost = all_cost / (float) flow_limit;
                }
            }
        }
        else{
            if (topoInfo.net_nodes[i].flow_ <= flow_limit){
                all_cost += topoInfo.net_nodes[i].flow_ * topoInfo.width_unit;
                topoInfo.net_nodes[i].unit_cost = all_cost / (float) topoInfo.net_nodes[i].flow_;
            }
            else {
                all_cost += flow_limit * topoInfo.width_unit;
                topoInfo.net_nodes[i].unit_cost = all_cost / (float) flow_limit;
            }
        }
    }
}

void middleCase2() {
    isScale = false;
    vector<ServerNode> init_server;
    for (int i = 0; i < NetNodeCount; i++) { // 直接链接的服务器
        if (!topoInfo.key_bit.test(i) && topoInfo.net_nodes[i].degree_ < 3) continue; // wh
        init_server.push_back(ServerNode(i, TopServerLevel));
    }
    unordered_map<int, Solution> init_solu;
    
    if (TopServerLevel + 2 > topoInfo.ServerLevelCount-1) {
        TopServerLevel = topoInfo.ServerLevelCount-1
    }else{
        TopServerLevel += 2;
    }
    
    for (int i = 0; i < 3; ++i){
        TopServerLevel -= i;
        TopServerFlow = topoInfo.serverLevel[TopServerLevel].outPut;
        TopServer = topoInfo.serverLevel[TopServerLevel];
        
        change();
        
        binaryinit(1, init_server);
        
        Solution s(g_servers, g_best_roads, totalcost_flow,g_servers_bit);
        s.toplevel = TopServerLevel;
        
        netnodeinit(s, 0);
        Solution ss(g_servers, g_best_roads, totalcost_flow,g_servers_bit);
        ss.toplevel = TopServerLevel;
        q.push(ss);
        
        TopServerLevel += i;
        
    }
//    for (int i = 0; i < 3; i++) {
//        TopServerLevel -= i;
//        netnodeinit(init_solu[i], 0);
//        TopServerLevel += i;
//        Solution s(g_servers, g_best_roads, totalcost_flow,g_servers_bit);
//        s.toplevel = TopServerLevel - i;
//        q.push(s);
//    }
    priority_queue<Solution> p = q;
    q = priority_queue<Solution>();
    int prevCost = MAX_NUMBER;
    int bestLevel = TopServerLevel;
    PRINT("find level \n");
    while (!TimeOutMin() && p.empty() == false) {
        Solution st = p.top();
        p.pop();
        g_servers = st.server;
        g_best_roads = st.roads;
        totalcost_flow = st.cost;
        g_servers_bit = st.bit;
        TopServerLevel = st.toplevel;
        PRINT("cur level %d \n", TopServerLevel);
        isScale = false;
        xjbs(0);
        xjbs_route(0);
        isScale = true;
        xjbdel();
        fireCore(g_servers, vector<ServerNode>());
        for (int i = 0; i < 8 && !TimeOutMin(); i++) {
            PRINT("-----------------step %d-----------------\n", i+1);
            if (TimeOutMin()) break;
            if (i % 2 == 0)
                xjbs(i);
            else
                xjbs_route(i);
            xjbdel();
            if (TimeOutMin()) break;
        }
        Solution s(g_servers, g_best_roads, totalcost_flow,g_servers_bit);
        s.toplevel = st.toplevel;
        q.push(s);
        if (TimeOutMin()) {
            PRINT("time out");
            break;
        }
    }
}
void highCase() {
    isScale = false;
    vector<ServerNode> init_server;
    for (int i = 0; i < NetNodeCount; i++) { // 直接链接的服务器
        if (!topoInfo.key_bit.test(i) && topoInfo.net_nodes[i].degree_ < 3) continue; // wh
        init_server.push_back(ServerNode(i, TopServerLevel));
    }
    unordered_map<int, Solution> init_solu;
    for (int i = 1; i < 4; ++i){
        binaryinit(i, init_server);
        Solution s(g_servers, g_best_roads, totalcost_flow,g_servers_bit);
        init_solu[i] = s;
    }
    for (int i = 1; i < 4; i++) {
        netnodeinit(init_solu[i], 0);
        q.push(Solution(g_servers, g_best_roads, totalcost_flow,g_servers_bit));
    }
    priority_queue<Solution> p = q;
    q = priority_queue<Solution>();
    int prevCost = MAX_NUMBER;
    while (!TimeOutMin() && p.empty()==false) {
        iter_size++;
        Solution st = p.top();
        p.pop();
        if (st.cost == prevCost) continue;
        prevCost = st.cost;
        g_servers = st.server;
        g_best_roads = st.roads;
        totalcost_flow = st.cost;
        g_servers_bit = st.bit;
        isScale = false;
        //xjbs_route(0);
        xjbs(0);
        isScale = true;
        xjbdel();
        fireCore(g_servers, vector<ServerNode>());
        for (int i = 0; i < 5 && !TimeOutMin(); i++) {
            PRINT("-----------------step %d-----------------\n", i+1);
            isScale = true;
            if (TimeOutMin()) break;
            xjbs(i);
        }
        //------------
        q.push(Solution(g_servers, g_best_roads, totalcost_flow,g_servers_bit));
        PRINT("*******totalcost:%d********\n",totalcost_flow);
    }
}



//  ------ find toplevel -----
// k遍历层数 len 表示某一层节点的数目， conumber表示是对对应的消费者的bfs，

void bfs_core(int k, bfsnode *b,int len, int conumber, int c, bitset<MaxNetNode> &usenode, queue<int>& node, priority_queue<bfsnode> bfs_server[MaxConsumNode]) {
    if (k<=0) {
        return ;
    }
    int n = 0;
    int root = 0;
    c = c+1;
    Graph p;
    for (int i = 0; i<len; ++i) {
        root = node.front();
        node.pop();
        if (!usenode.test(root)) {
            b->node.push_back(root);
            usenode.set(root);
        }
        p = topoInfo.graph[root];
        while(p) {
            node.push(p->endPoint);
            n = n+1;
            p = p->next;
        }
    }
    k--;
    b->pathcost = c;
    bfs_server[conumber].push(*b);
    b->node.clear();
    bfs_core(k,b,n, conumber,  c, usenode, node,bfs_server);
}


bool comp_d_f_c(const int &a, const int &b){ // 给网络节点排序，，先度 后输出能力 再成本
    NetNode a_ = topoInfo.net_nodes[a];
    NetNode b_ = topoInfo.net_nodes[b];
    
    if (a_.degree_ == b_.degree_) {
        if (a_.flow_ ==  b_.flow_) {
            return a_.cost < b_.cost;
        }
        return a_.flow_ > b_.flow_;
    }
    return a_.degree_ > b_.degree_;
}
/**
 统计每个消费节点距离为k 的节点的频次
 获取 count 个
 
 @param k 距离
 @param count 个数
 @return 结果
 */
vector<ServerNode> mybfs(int k, int count) {
    bfsnode b;
    priority_queue<bfsnode> bfs_server[MaxConsumNode];
    bitset<MaxNetNode> usenode; // 记录某个节点是否已被遍历
    queue<int> node; //即将被遍历的点
    int c = 0; // 记录层数pathcost
    
    for (int i = 0; i<topoInfo.ConsumNodeCount; ++i) {
        while(!node.empty()) {
            node.pop();
        }
        c = 0;

        usenode.reset();
        node.push(topoInfo.consumNode[i].netNode);
        bfs_core(k, &b, 1, i, c, usenode, node,bfs_server);
    }
    // 获取度大的。
    unordered_set<int> aa;
    for (int i = 0; i < ConsumNodeCount; i++) {
        bfsnode a;
        a = bfs_server[i].top();
        for (int i = 0; i < a.node.size(); i++) {
            aa.insert(a.node[i]);
        }
    }
    
    vector<int> addup;
    addup.assign(aa.begin(), aa.end());
    
    sort(addup.begin(), addup.end(), comp_d_f_c);
    vector<ServerNode> getnode;
    
    for (int i = 0; i < MIN(count, aa.size()); i++) {
        getnode.push_back(ServerNode(addup[i], topoInfo.net_nodes[addup[i]].maxServerLevel));
    }
    //PRINT("\n");
    return getnode;
}
bool c(const ServerNode &a, const ServerNode &b)
{
    if (topoInfo.net_nodes[a.netnode].degree_ == topoInfo.net_nodes[b.netnode].degree_){
        return topoInfo.net_nodes[a.netnode].flow_ < topoInfo.net_nodes[b.netnode].flow_;
    }
    return topoInfo.net_nodes[a.netnode].degree_ < topoInfo.net_nodes[b.netnode].degree_;
}

void find_topleve_lbinary() {
    totalcost_flow = MAX_NUMBER;
    vector<ServerNode> server,tmp,tmp2,back;
    for (int i = 0; i < NetNodeCount; i++) { // 直接链接的服务器
        if (!topoInfo.key_bit.test(i) && topoInfo.net_nodes[i].degree_ < 3) continue; // wh
        
        server.push_back(ServerNode(i, TopServerLevel));
    }
    
    sort(server.begin(),server.end(),c);
    
    int num = 1;

    while (true) {
        tmp.clear();
        for (int i = num; i < server.size(); i++) {
            tmp.push_back(server[i]);
        }
        tmp2 = tmp;
        for (int i = 0 ;i < back.size(); i++)
            tmp2.push_back(back[i]);
        
        if (findMaxFlow(tmp2)) {
            server = tmp;
            num *= 2;
        } else {
            if (num != 1) {
                num = 1;
            } else {
                back.push_back(server.front());
            }
        }
        if (num >= server.size()) {
            break;
        }
    }
    findRoadsProgrammeWithServerNodebymaxflow(tmp2);
}
void find_toplevel_init(){
    int size;
    totalcost_flow = MAX_NUMBER;
    bitset<MaxNetNode> bit;
    vector<ServerNode> server;
    
    find_topleve_lbinary();
    server = g_servers;
    bit = g_servers_bit;

    //  添加 --> 删除
    for (int i = 0; i < 3; i++) {
        int n = 50 - i*22;
        // 添加
        vector<ServerNode> choosenode = mybfs(3,50);

        PRINT("choosenode: %d\n", choosenode.size());
        for (int i = 0; i < MIN(n, choosenode.size()); i++) {
            if (!bit.test(choosenode[i].netnode)) { // 服务器没 加
                choosenode[i].serverLevel = TopServerLevel;
                server.push_back(choosenode[i]);
            }
        }

        //for
        size = (int)server.size();
        for (int i = 0;i < size*1; i++) {
            int index = randomi(0,(int)server.size() - 1);
            ServerNode temp = server[index];
            server.erase(server.begin()+index);
           
            if (!findMaxFlow(server)) {// 删除失败
                server.push_back(temp);
            }
        }

        findRoadsProgrammeWithServerNodebymaxflow(server);
        PRINT("serversize:%d\ttotalcost:%d\n", (int)server.size(), totalcost_flow);
        
    }//for
    
}

//你要完成的功能总入口
void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    TimeOut();
    TimeOutMin();
    initial(topo, line_num, topoInfo);
    init();
    
    int n;
    if (NetNodeCount > 1000) {
        n = 6;
    }else{
        n = 7;
    }
    
    isScale = false;
    PRINT("Begin Init:%d\n", RunTime());
    isfindtoplevel = true;
    priority_queue<Solution> initq;
    for (int i = 0; i < n; i++) {
        if (TopServerLevel - i < 0) {
            break;
        }
        
        TopServerLevel -= i;
        PRINT("%d ---\n ", TopServerLevel);
        TopServerFlow = topoInfo.serverLevel[TopServerLevel].outPut;
        TopServer = topoInfo.serverLevel[TopServerLevel];
        
        find_toplevel_init();
        // test
        //isScale = true;
        //fireCore(g_servers, vector<ServerNode>());
        //isScale = false;
        initq.push(Solution(g_servers, g_best_roads, totalcost_flow,TopServerLevel));
        
        TopServerLevel += i;
        TopServerFlow = topoInfo.serverLevel[TopServerLevel].outPut;
        TopServer = topoInfo.serverLevel[TopServerLevel];
    }
    
    Solution s = initq.top();
    
    // 设置最高档次
    PRINT("TopLevel: %d\n", s.toplevel);

    TopServerLevel = s.toplevel;
    // if (NetNodeCount > 800) 
    //TopServerLevel += 1;
    TopServerFlow = topoInfo.serverLevel[TopServerLevel].outPut;
    TopServer = topoInfo.serverLevel[TopServerLevel];
    PRINT("End Init:%d\n", RunTime());

    isfindtoplevel = false;

    /* 计算单位成本和服务器最高档次 */
    int all_cost, flow_limit;
    for (int i = 0; i < NetNodeCount; i++) {
        all_cost = 0;
        flow_limit = 0;

        int output = topoInfo.net_nodes[i].flow_;
        if(output >= TopServerFlow){
            topoInfo.net_nodes[i].maxServerLevel = TopServerLevel;
            all_cost += TopServer.cost;
            flow_limit = TopServerFlow;
        } else {
            ServerLevel s = checkServerLevel(output);
            topoInfo.net_nodes[i].maxServerLevel = s.levelid;
            all_cost += s.cost;
            flow_limit = s.outPut;
        }

        all_cost += topoInfo.net_nodes[i].cost;

        if (topoInfo.key_bit.test(i)) {
            int cid = topoInfo.g_net2consume[i];
            if (topoInfo.net_nodes[i].flow_ <= flow_limit){
                all_cost += (topoInfo.net_nodes[i].flow_ - topoInfo.consumNode[cid].need) * topoInfo.width_unit;
                topoInfo.net_nodes[i].unit_cost = all_cost / (float) topoInfo.net_nodes[i].flow_;
            }
            else{
                if(topoInfo.consumNode[cid].need < flow_limit){
                    all_cost += (flow_limit - topoInfo.consumNode[cid].need) * topoInfo.width_unit;
                    topoInfo.net_nodes[i].unit_cost = all_cost / (float) flow_limit;
                }
                else{
                     topoInfo.net_nodes[i].unit_cost = all_cost / (float) flow_limit;
                }
            }
        }
        else{
            if (topoInfo.net_nodes[i].flow_ <= flow_limit){
                all_cost += topoInfo.net_nodes[i].flow_ * topoInfo.width_unit;
                topoInfo.net_nodes[i].unit_cost = all_cost / (float) topoInfo.net_nodes[i].flow_;
            }
            else {
                all_cost += flow_limit * topoInfo.width_unit;
                topoInfo.net_nodes[i].unit_cost = all_cost / (float) flow_limit;
            }
        }
    }
    
    /* 构造初始解 */


    // 核心计算

    if (NetNodeCount < 800){
        TIME_OUTMIN = 80;
        netnodeinit_ccc = 0.0019;
        g_slave = 1.00031;
        //middleCase();
        middleCase2();
    }else{
        TIME_OUTMIN = 80;
        netnodeinit_ccc = 0.0017;
        g_slave = 1.00010;
        highCase();
    }
    
    isScale = true;
    
    Solution st = q.top();
    g_servers = st.server;
    g_best_roads = st.roads;
    totalcost_flow = st.cost;
    g_servers_bit = st.bit;
    
    while (!TimeOut()) {
        fireCore(g_servers, vector<ServerNode>());
    }
    sort(g_servers.begin(),g_servers.end(),servercmp);
    
    PRINT("sucfind: %d\n", suc_find);
    PRINT("iter time %d \n",iter_size);
    PRINT("server size %d \n",(int)g_servers.size());
    writeresult(filename); 

    PRINT("Finish:%d\n", RunTime());
}
