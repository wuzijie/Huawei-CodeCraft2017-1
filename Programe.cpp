//
//  Programe.cpp
//  CDN复赛
//
//  Created by Daniel on 2017/4/9.
//  Copyright © 2017年 Daniel. All rights reserved.
//

#include "cdn.h"


/**
 输出用最大流获取的结果
 
 @param result 存放输出字符串
 @param g_servers 服务器
 @param g_best_roads 使用到的链路集合
 @param NetNodeCount 网络节点个数
 @param g_net2consume 网络节点id转消费节点id
 */
void writeMaxFlowPrograme(string& result, vector<ServerNode>& g_servers, vector<Edge>& g_best_roads, int NetNodeCount, vector<int>& g_net2consume) {
    
    // 存储路径
    unordered_map<int,vector<Edge> > road_map;
    unordered_map<int,Edge> other_road_map;
    unordered_map<int, int> serverlevel;
    
    bitset<MaxNetNode> b;
    for (int i = 0; i < g_servers.size(); i++) {
        b.set(g_servers[i].netnode);
    }
    for (int i = 0; i < g_best_roads.size(); i++) {
        if (b[g_best_roads[i].end] == 1 && g_best_roads[i].start == NetNodeCount) {
            ServerLevel s = checkServerLevel(g_best_roads[i].flow);
            serverlevel[g_best_roads[i].end] = s.levelid;
        }
    }
    
    
    for (int i = 0; i < g_best_roads.size(); i++) {
        Edge e = g_best_roads[i];
        if (road_map.find(e.start)==road_map.end())  // 寻找路径用
            road_map[e.start] = vector<Edge>();
        road_map[e.start].push_back(e);
        
        other_road_map[e.start*10000+e.end] = e; // 分配路径用
    }
    queue<Flow> q;
    Flow f;
    f.cur_node = NetNodeCount;
    f.need = 999999;
    f.roads.push_back(NetNodeCount);
    q.push(f); // 超级源
    vector<string> answers;     // 寻找出来路径的 string 后面已加 消费者。
    vector<vector<int>> roads;  // 寻找出来的路径  最后一个是超级汇点  下标与 answers。对应
    while (!q.empty()) {
        //        printf("%d\n",(int)q.size());
        Flow f = q.front();
        q.pop();
        // 如果没有后路，当前节点是消费节点对应的网络节点
        if (road_map[f.cur_node].size() == 0) { // 到超汇点
            
            answers.push_back(f.road);  //记录 路径 string
            roads.push_back(f.roads);   //记录 路径 用于分配带宽
        }
        //PRINT("size : %d\n",road_map[f.cur_node].size());
        vector<Edge> v = road_map[f.cur_node];
        for (int i = 0; i < v.size(); i++) {
            
            //PRINT("%d to %d\n",f.cur_node,v[i].t);
            Flow tmp = f;
            char c[20];
            
            if (v[i].end == NetNodeCount+1) { // 最后一个 输入 消费者。
                sprintf(c, "%d ",g_net2consume[v[i].start]);
                tmp.roads.push_back(v[i].end); // 最后一个加入超汇点
                
            }else{ // 输入网络节点
                sprintf(c, "%d ",v[i].end);
                tmp.roads.push_back(v[i].end);
            }
            tmp.cur_node = v[i].end;
            tmp.road += c;
            q.push(tmp);
        }
    }
    
    string road_string = "";
    int road_num = 0;
    for (int i = 0; i < answers.size(); i++) {  /// 在这里确定 路径的 宽带。。 路径上的链路有些为0的时候就删除这条路
        // check 每条路。
        vector<int> r = roads[i];
        
        int minWidth = 99999999;
        for (int j = 0; j < r.size()-1; j++) {
            int start   = r[j];
            int end     = r[j+1];
            if (other_road_map[start*10000+end].flow < minWidth) { // 找最小
                minWidth = other_road_map[start*10000+end].flow;
            }
        }
        if (minWidth != 0) { // 可以分配 带宽
            string s = answers[i] + to_string(minWidth) + " " + to_string(serverlevel[r[1]]) + "\n";    // 末尾加上 带宽
            road_string += s;
            road_num++;
            
            for (int j = 0; j < r.size()-1; j++) { // 减去。使用过的
                int start   = r[j];
                int end     = r[j+1];
                other_road_map[start*10000+end].flow -= minWidth;
            }
        }
    }
    char c[20];
    
    sprintf(c, "%d\n\n", road_num);
    result = c + road_string;
    
    result = result.substr(0,result.size()-1); //去掉最后一个 \n
    
    return;
}
