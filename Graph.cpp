//
//  Graph.cpp
//  CDN复赛
//
//  Created by Daniel on 2017/4/9.
//  Copyright © 2017年 Daniel. All rights reserved.
//

#include "cdn.h"


/**
 链路 反向
 
 @param link1 原链路指针
 @return 新链路指针
 */
NetLink * RNetLink(NetLink *link1){
    NetLink *link2 = (NetLink *)malloc(sizeof(NetLink));
    link2->startPoint   = link1->endPoint;
    link2->endPoint     = link1->startPoint;
    link2->maxWidth     = link1->maxWidth;
    link2->remainWidth  = link1->remainWidth;
    link2->cost         = link1->cost;
    link2->next         = NULL;
    return link2;
}
// 设置边
void Mcmf::set_arc(int tail_id, int head_id, int capacity,int cost){
    arc_first[tail_id + 1] ++;
    arc_first[head_id + 1] ++;
    i_node = nodes + tail_id;
    j_node = nodes + head_id;
    
    int pos_sister = (int)pos_current + 1;
    McmfEdge* arc_sister = arc_current + 1;
    
    arc_tail[pos_current]   = tail_id;
    arc_tail[pos_sister] = head_id;
    arc_current->head = j_node ;
    arc_current->tail = i_node;
    arc_current->rez_capacity = capacity ;
    cap[pos_current] = capacity;
    arc_current->cost = cost ;
    arc_current->sister = arc_sister ;
    
    arc_sister->head = i_node;
    arc_sister->tail = j_node;
    arc_sister->rez_capacity =0 ;
    cap[pos_sister] = 0;
    arc_sister->cost = -cost ;
    arc_sister->sister = arc_current ;
    
    if (cost > max_cost) {
        max_cost = cost;
    }
    
    arc_current += 2;
    pos_current += 2;
}



/**
 图的初始化
 
 @param topo 图文件
 @param line_num 行数
 @param topoInfo 图基本信息
 @return 0
 */
int initial(char * topo[MAX_EDGE_NUM], int line_num, TopoInfo& topoInfo){
    int number;

    int all_cost = 0;
    int all_width = 0;

    char *c; // 获取字符串
    // 1.第一行：网络节点数量 网络链路数量 消费节点数量
    c = topo[0];
    sscanf(c, "%d%d%d", &topoInfo.NetNodeCount, &topoInfo.LinkCount, &topoInfo.ConsumNodeCount);
    
    // 初始化网络节点标
    topoInfo.net_nodes.resize(topoInfo.NetNodeCount, NetNode());
    topoInfo.consumNode.resize(topoInfo.ConsumNodeCount, ConsumNode());
    topoInfo.g_net2consume.resize(topoInfo.NetNodeCount, -1);
    
    
    // 2. 服务器档次 成本
    for (number = 2; number < line_num; number++) {
        c = topo[number];
        if (*c == '\n' || *c == '\r' || *c == '\0') {
            break;
        }
        ServerLevel s;
        sscanf(c, "%d%d%d", &s.levelid,&s.outPut,&s.cost);
        topoInfo.serverLevel.push_back(s);
    }
    topoInfo.ServerLevelCount = (int)topoInfo.serverLevel.size();
    // 3. 网络节点成本
    for (number = number + 1; number < line_num; number++) {
        c = topo[number];
        if (*c == '\n' || *c == '\r' || *c == '\0') {
            break;
        }
        int id,cost;
        sscanf(c, "%d%d", &id, &cost);
        topoInfo.net_nodes[id].cost = cost;
    }
    
    // 4. 链路 图
    for (int i = 0; i < topoInfo.NetNodeCount; i++) {
        Graph a = NULL;
        topoInfo.graph.push_back(a);
    }
    
    int start_p, end_p, width, cost;
    NetLink *link;
    for (number = number + 1; number < line_num; number++) {
        c = topo[number];
        sscanf(c, "%d%d%d%d", &start_p, &end_p, &width, &cost);
        if (*c == '\n' || *c == '\r' || *c == '\0') {
            break;
        }
        //更新网络节点表
        topoInfo.net_nodes[start_p].flow_ += width;
        topoInfo.net_nodes[end_p].flow_ += width;
        
        topoInfo.net_nodes[start_p].degree_++;
        topoInfo.net_nodes[end_p].degree_++;
        
        // 最大流需要用边
        topoInfo.g_edges_maxflow.push_back(NetLink(start_p,end_p,width,cost));
        
        link = (NetLink *)malloc(sizeof(NetLink));
        link->startPoint = start_p;
        link->endPoint = end_p;
        link->maxWidth = width;
        link->remainWidth = width;
        link->cost = cost;
        link->next = NULL;
        
        link->next = topoInfo.graph[start_p];
        topoInfo.graph[start_p] = link;
        // 反向 存另外一条
        link = RNetLink(link);
        link->next = topoInfo.graph[end_p];
        topoInfo.graph[end_p] = link;
        
        all_cost += cost * width;
        all_width += width;
    }// for
    
    topoInfo.all_demend = 0;
    // 4. 消费节点 储存
    int netnode, need, j;
    for (number = number+1; number < line_num; number++) {
        c = topo[number];
        sscanf(c, "%d%d%d", &j, &netnode, &need);
        topoInfo.consumNode[j].netNode  = netnode;
        topoInfo.consumNode[j].need     = need;
        topoInfo.g_net2consume[netnode] = j;
        topoInfo.all_demend += need;
        
        topoInfo.net_nodes[netnode].flow_ += need;
        topoInfo.key_bit.set(netnode);
    }//for
    
    topoInfo.width_unit = ceil(all_cost / (float) all_width * 1.7);

    return 1;
}






int REDUCED_COST(McmfNode *i, McmfNode *j, McmfEdge *a  ) {
    return i->price + a->cost - j->price;
}

void Mcmf::allocate_arrays(){
    //给边和节点分配内存
    nodes = (McmfNode*) calloc ( n+2,   sizeof(McmfNode) );
    arcs = (McmfEdge*)  calloc ( 2*m+1, sizeof(McmfEdge) );
    cap = (int*) calloc ( 2*m,   sizeof(int) );
    
    arc_tail = (int*) calloc ( 2*m,   sizeof(int) );
    arc_first = (int*) calloc ( n+2,   sizeof(int) );
    
    for ( McmfNode *in = nodes; in <= nodes + n; in ++ ) {
        in->excess = 0;
    }
    if ( nodes == NULL || arcs == NULL || arc_first == NULL || arc_tail == NULL) {
        //		PRINT("Error:  Memory allocation problem inside algo\n");
    }
    
    pos_current = 0;
    arc_current = arcs;
    
    max_cost = 0;
    total_p = total_n = 0;
}
// 扫面节点
void Mcmf::up_node_scan( McmfNode *i){
    int j_new_rank;
    int i_rank = i->rank;
    
    // 扫描边
    McmfEdge *a = i->first;
    McmfEdge *a_stop = (i + 1)->suspended;
    while (a != a_stop) {
        McmfEdge *ra = a->sister;
        if (ra->rez_capacity > 0) {
            McmfNode *j = a->head;
            int j_rank = j->rank;
            if ( j_rank > i_rank ) {
                int rc = REDUCED_COST( j, i, ra );
                if (rc < 0) {
                    j_new_rank = i_rank;
                } else {
                    int dr = rc / epsilon;
                    j_new_rank = ( dr < linf ) ? i_rank + dr + 1 : linf;
                }
                if ( j_rank > j_new_rank ) {
                    j->rank = j_new_rank;
                    j->current = ra;
                    if ( j_rank < linf ) {
                        McmfBucket *b_old = buckets + j_rank;
                        remove_from_bucket( j, b_old );
                    }
                    McmfBucket *b_new = buckets + j_new_rank;
                    insert_to_bucket( j, b_new );
                }
            }
        }
        a++;
    }
    
    i->price -=  i_rank * epsilon;
    i->rank = -1;
}

// 删除分配的数组
void Mcmf::deallocate_arrays(){
    if ( arcs) free ( arcs );
    if ( dnode) delete dnode;
    if ( cap) free ( cap );
    if ( buckets) free ( buckets );
    
    if ( nodes) {
        nodes = nodes - node_min;
        free ( nodes );
    }
}
// 更新约束
int Mcmf::update_epsilon(){
    if (epsilon <= 1){
        return 1;
    }
    epsilon = ceil((double)epsilon/f_scale);
    cut_off = cut_off_factor * epsilon;
    cut_on = cut_off * 0.8;
    return 0;
}


// 设置提供需求的节点
void Mcmf::set_supply_demand_of_node( int head_id, int excess){
    (nodes + head_id)->excess = excess;
    if ( excess > 0) total_p += excess;
    if ( excess < 0) total_n -= excess;
}

// 预处理
void Mcmf::pre_processing(){
    int last, arc_num, arc_new_num;;
    int tail_id;
    McmfNode *head_p, *tail_p;
    McmfEdge *arc_new, *arc_tmp;
    int up_bound;
    int cost;
    
    (nodes + node_min)->first = arcs;
    
    for (int  i = node_min + 1; i <= node_max + 1; i ++ ) {
        arc_first[i] += arc_first[i-1];
        ( nodes + i )->first = arcs + arc_first[i];
    }
    
    for (int  i = node_min; i < node_max; i ++ ) {
        last = (int)((( nodes + i + 1 )->first) - arcs);
        for ( arc_num = arc_first[i]; arc_num < last; arc_num ++ ) {
            tail_id = arc_tail[arc_num];
            while ( tail_id != i ) {
                arc_new_num = arc_first[tail_id];
                arc_current = arcs + arc_num;
                arc_new = arcs + arc_new_num;
                
                head_p = arc_new->head;
                arc_new->head = arc_current->head ;
                arc_current->head = head_p ;
                
                tail_p = arc_new->tail;
                arc_new->tail =  arc_current->tail;
                arc_current->tail = tail_p;
                
                up_bound          = cap[arc_new_num];
                cap[arc_new_num] = cap[arc_num];
                cap[arc_num]     = up_bound;
                
                up_bound = arc_new->rez_capacity;
                arc_new->rez_capacity = arc_current->rez_capacity ;
                arc_current->rez_capacity = up_bound ;
                
                cost = arc_new->cost;
                arc_new->cost = arc_current->cost ;
                arc_current->cost = cost ;
                
                if ( arc_new != arc_current->sister ) {
                    arc_tmp = arc_new->sister;
                    arc_new->sister = arc_current->sister ;
                    arc_current->sister = arc_tmp ;
                    
                    arc_current->sister->sister = arc_current ;
                    arc_new->sister->sister = arc_new ;
                }
                
                arc_tail[arc_num] = arc_tail[arc_new_num];
                arc_tail[arc_new_num] = tail_id;
                
                arc_first[tail_id] ++ ;
                
                tail_id = arc_tail[arc_num];
            }
        }
    }
    
    n = node_max - node_min + 1;
    nodes = nodes + node_min;
    
    free ( arc_first );
    free ( arc_tail );
}

// 初始化
void Mcmf::algo_initialize(){
    
    f_scale = (long) 12;
    sentinel_node = nodes + n;
    sentinel_arc = arcs + m;
    
    for (McmfNode *i = nodes; i != sentinel_node; i++) {
        i->price = 0;
        i->suspended = i->first;
        i->q_next = sentinel_node;
    }
    sentinel_node->first = sentinel_arc;
    sentinel_node->suspended =  sentinel_arc;
    
    for (McmfNode *i = nodes; i != sentinel_node; i++) {
        McmfEdge *a = i->first;
        McmfEdge *a_stop = (i + 1)->suspended;
        while (a != a_stop){
            if ( a->cost < 0) {
                long df = a->rez_capacity;
                if (df > 0) {
                    increase_flow( i, a->head, a, df);
                }
            }
            a++;
        }
    }
    dn = n + 1;
    for (McmfEdge *a = arcs; a != sentinel_arc; a ++ ) {
        a->cost *=  dn;
    }
    mmc = max_cost * dn;
    linf = (dn * ceil(f_scale) + 2);
    buckets = (McmfBucket*) calloc ( linf, sizeof(McmfBucket));
    if ( buckets == NULL ) {
        this->is_error = true;
        //		PRINT("ALLOCATION_FAULT");
        return;
    }
    l_bucket = buckets + linf;
    dnode = new McmfNode;
    for (McmfBucket *b = buckets; b != l_bucket; b++) {
        b->p_first = dnode;
    }
    epsilon = mmc;
    if ( epsilon < 1) {
        epsilon = 1;
    }
    price_min = -MAX_NUMBER;
    cut_off_factor = 1.5 * pow( (double)n, 0.44);
    cut_off_factor = MAX( cut_off_factor, 12);
    n_ref = 0;
    flag_price = 0;
    excq_first = NULL;
}



// 更新price
void Mcmf::price_update(){
    for (McmfNode *i = nodes; i != sentinel_node; i ++ ) {
        if ( i->excess < 0 ) {
            insert_to_bucket( i, buckets );
            i->rank = 0;
        } else {
            i->rank = linf;
        }
    }
    
    long remain = total_excess;
    if (remain < 0.5)
        return;
    
    McmfBucket *b;
    for (b = buckets; b != l_bucket; b++) {
        
        while ( b->p_first != dnode ) {
            McmfNode *i = b->p_first;
            b->p_first = i->b_next;
            
            up_node_scan( i );
            
            if ( i ->excess > 0 ) {
                remain -= ( i->excess);
                if ( remain <= 0 ) break;
            }
        }
        if ( remain <= 0 ) break;
    }
    
    if (remain > 0.5)
        flag_updt = 1;
    
    long dp = (b - buckets) * epsilon;
    for ( McmfNode *i = nodes; i != sentinel_node; i ++ ) {
        if ( i->rank >= 0 ) {
            if ( i->rank < linf ) {
                remove_from_bucket( i, ( buckets + i->rank) );
            }
            if ( i->price > price_min ) {
                i->price -=  dp;
            }
        }
    }
}

// 重标签
int Mcmf::relabel( McmfNode *i){
    int p_max = price_min;
    int i_price = i->price;
    
    McmfEdge *a_max = NULL;
    
    McmfEdge *a = i->current + 1;
    McmfEdge *a_stop = (i + 1)->suspended;
    while (a != a_stop) {
        int dp = (a->head->price - a->cost);
        if ( a->rez_capacity > 0 && ( dp > p_max ) ) {
            if ( i_price < dp ) {
                i->current = a;
                return 1;
            }
            p_max = dp;
            a_max = a;
        }
        a++;
    }
    
    for ( a = i->first, a_stop = i->current + 1; a != a_stop; a ++ ) {
        int dp = (a->head->price - a->cost);
        if ( a->rez_capacity > 0 && ( dp > p_max ) ) {
            if ( i_price < dp ) {
                i->current = a;
                return 1;
            }
            p_max = dp;
            a_max = a;
        }
    }
    
    if ( p_max != price_min ) {
        i->price =  p_max - epsilon;
        i->current = a_max;
    }
    else {
        if ( i->suspended == i->first ) {
            if ( i->excess == 0 ) {
                i->price =  price_min;
            } else {
                if ( n_ref == 1 ) {
                    this -> is_error = true;
                    //					PRINT("UNFEASIBLE");
                    return 0;
                } else {
                    this -> is_error = true;
                    //					PRINT("PRICE_OFL");
                    return 0;
                }
            }
        } else {
            flag_price = 1;
        }
    }
    
    n_rel ++;
    return 0;
}

void Mcmf::discharge( McmfNode *i){
    McmfEdge *a = i->current;
    McmfNode *j = a->head;
    
    if (!((a->rez_capacity > 0) && (i->price + a->cost < j->price))) {
        relabel( i );
        if (this -> is_error) return;
        a = i->current;
        j = a->head;
    }
    
    while (1) {
        long j_exc = j->excess;
        if ( j_exc >= 0 ) {
            long df = MIN( i->excess, a->rez_capacity );
            if (j_exc == 0){
                n_src++;
            }
            increase_flow( i, j, a, df );
            if (j ->q_next == sentinel_node) {
                insert_to_excess_q( j );
            }
        }
        else {
            long df = MIN( i->excess, a->rez_capacity );
            increase_flow( i, j, a, df );
            if ( j->excess >= 0 ) {
                if ( j->excess > 0 ) {
                    n_src ++;
                    relabel( j );
                    if (this -> is_error) return;
                    insert_to_excess_q( j );
                }
                total_excess += j_exc;
            }
            else {
                total_excess -= df;
            }
        }
        
        if ( i->excess <= 0) n_src --;
        if ( i->excess <= 0 || flag_price ) break;
        
        relabel( i );
        if (this -> is_error) return;
        a = i->current;
        j = a->head;
    }
    i->current = a;
}

int Mcmf::price_in(){
    int bad_found = 0;
    int n_in_bad = 0;
    
    int n = 1;
    while (n == 1) {
        n = 0;
        for (McmfNode *i = nodes; i != sentinel_node; i ++ ) {
            McmfEdge *a = i->first - 1;
            McmfEdge *a_stop = i->suspended - 1;
            while (a != a_stop) {
                long rc = REDUCED_COST( i, a->head, a );
                if ( ( rc < 0) && ( a->rez_capacity > 0) ) { // 坏的用例;
                    if ( bad_found == 0 ) {
                        bad_found = 1;
                        update_cut_off();
                        n = 1;
                        break;
                    }
                    increase_flow( i, a->head, a, a->rez_capacity );
                    
                    McmfEdge *ra = a->sister;
                    McmfNode *j  = a->head;
                    
                    i->first--;
                    exchange( a, i->first );
                    if ( ra < j->first ) {
                        j->first--;
                        exchange(ra, j->first);
                    }
                    
                    n_in_bad ++;
                }
                else {
                    if ( ( rc < cut_on ) && ( rc > -cut_on ) ) {
                        i->first--;
                        exchange( a, i->first );
                    }
                }
                a--;
            }//while
            if (n == 1) {
                break;
            }
        }//for
    }//while
    
    if (n_in_bad != 0) {
        n_bad_pricein ++;
        total_excess = 0;
        n_src = 0;
        reset_excess_q();
        
        for (McmfNode *i = nodes; i != sentinel_node; i ++ ) {
            i->current = i->first;
            long i_exc = i->excess;
            if (i_exc > 0) {
                
                total_excess += i_exc;
                n_src ++;
                insert_to_excess_q( i );
            }
        }
        
    }
    if ( time_for_price_in == 4){
        time_for_price_in = 6;
    }
    if ( time_for_price_in == 2){
        time_for_price_in = 4;
    }
    return n_in_bad;
}


void Mcmf::refine(){
    int pr_in_int;
    
    n_ref ++;
    n_rel = 0;
    pr_in_int = 0;
    // 初始化
    total_excess = 0;
    n_src = 0;
    reset_excess_q();
    time_for_price_in = 2;
    
    for ( McmfNode *i = nodes; i != sentinel_node; i ++ ) {
        i->current = i->first;
        long i_exc = i->excess;
        if ( i_exc > 0 ) {
            total_excess += i_exc;
            n_src++;
            insert_to_excess_q(i);
        }
    }
    if (total_excess <= 0)
        return;
    
    while (1) {
        if (excq_first == NULL) {
            if ( n_ref > 1 ) {
                pr_in_int = 0;
                price_in();
            }
            if (excq_first == NULL)
                break;
        }
        McmfNode *i = excq_first;
        excq_first = i->q_next;
        i->q_next = sentinel_node;
        
        if ( i->excess > 0 ) {
            discharge( i );
            if (this -> is_error)
                return;
            if ( (n_rel > n * 0.4 + n_src * 30) || flag_price ) {
                if ( i->excess > 0 ) {
                    insert_to_excess_q( i );
                }
                if ( flag_price && ( n_ref > 1 ) ) {
                    pr_in_int = 0;
                    price_in();
                    flag_price = 0;
                }
                price_update();
                while (flag_updt) {
                    if (n_ref == 1) {
                        this-> is_error = true;
                        //						PRINT("UNFEASIBLE");
                        return;
                    }else{
                        flag_updt = 0;
                        update_cut_off();
                        n_bad_relabel ++;
                        pr_in_int = 0;
                        price_in();
                        price_update();
                    }
                }
                n_rel = 0;
                if (n_ref > 1 && (pr_in_int ++ > time_for_price_in)) {
                    pr_in_int = 0;
                    price_in();
                }
            }
        }
    }
    
    return;
}

int Mcmf::price_refine(){
    int cc = 1;
    int snc = 0;
    
    while (1) {
        int nnc = 0;
        for (McmfNode *i = nodes; i != sentinel_node; i ++ ) {
            i->rank = 0;
            i->inp = 0;
            i->current = i->first;
        }
        reset_excess_q();
        for (McmfNode *i = nodes; i != sentinel_node; i ++ ) {
            if (i->inp == 2)
                continue;
            i->b_next = NULL;
            while (1) {
                i->inp = 1;
                McmfEdge *a = i->current;
                McmfEdge *a_stop = (i + 1)->suspended;
                while (a != a_stop) {
                    if (a->rez_capacity > 0) {
                        McmfNode *j = a->head;
                        if (REDUCED_COST ( i, j, a ) < 0) {
                            if (j->inp == 0) {
                                i->current = a;
                                j->b_next = i;
                                i = j;
                                a = j->current;
                                a_stop = (j+1)->suspended;
                                break;
                            }
                            if (j->inp == 1) {
                                cc = 0;
                                nnc ++;
                                i->current = a;
                                McmfNode *is, *ir;
                                is = ir = i;
                                int df = MAX_NUMBER;
                                while (1) {
                                    McmfEdge *ar = ir->current;
                                    if (ar->rez_capacity <= df) {
                                        df = ar->rez_capacity;
                                        is = ir;
                                    }
                                    if (ir == j)
                                        break;
                                    ir = ir->b_next;
                                }//while
                                ir = i;
                                while (1) {
                                    McmfEdge *ar = ir->current;
                                    increase_flow( ir, ar->head, ar, df);
                                    if ( ir == j ) break;
                                    ir = ir->b_next;
                                }
                                if (is != i) {
                                    for ( ir = i; ir != is; ir = ir->b_next ) {
                                        ir->inp = 0;
                                    }
                                    i = is;
                                    a = is->current + 1;
                                    a_stop = (is+1)->suspended;
                                    break;
                                }//if
                            }//if
                        }//if
                    }//if
                    a++;
                }//while
                if (a == a_stop) {
                    i->inp = 2;
                    McmfNode *j = i->b_next;
                    i->q_next = excq_first;
                    excq_first = i;
                    if (j == NULL)
                        break;
                    i = j;
                    i->current++;
                }//if
            }//while
        }//for
        
        snc += nnc;
        if (snc < 0) {
            cc = 1;
        }
        if (cc == 0){
            break;
        }
        int bmax = 0;
        
        while (excq_first != NULL) {
            McmfNode *i = excq_first;
            excq_first = i->q_next;
            i->q_next = sentinel_node;
            int i_rank = i->rank;
            
            McmfEdge *a = i->first;
            McmfEdge *a_stop = (i + 1)->suspended;
            while (a != a_stop) {
                if ( a->rez_capacity > 0 ) {
                    McmfNode *j  = a->head;
                    int rc = REDUCED_COST( i, j, a );
                    
                    if ( rc < 0 ) {
                        int dr = (( - rc - 0.5 ) / epsilon);
                        int j_rank = dr + i_rank;
                        if (j_rank < linf ) {
                            if ( j_rank > j->rank )
                                j->rank = j_rank;
                        }//if
                    }//if
                }//if
                a++;
            }//while
            if (i_rank > 0) {
                if ( i_rank > bmax ) bmax = i_rank;
                McmfBucket *b = buckets + i_rank;
                insert_to_bucket( i, b );
            }
        }
        if (bmax == 0){
            break;
        }
        for (McmfBucket * b = buckets + bmax; b != buckets; b-- ) {
            int i_rank = (int)(b - buckets);
            int dp = i_rank * epsilon;
            
            while ( b->p_first != dnode ) {
                McmfNode *i = b->p_first;
                b->p_first = i->b_next;
                
                McmfEdge *a = i->first;
                McmfEdge *a_stop = (i + 1)->suspended;
                while (a != a_stop) {
                    if (a->rez_capacity > 0) {
                        McmfNode *j = a->head;
                        int j_rank = j->rank;
                        if (j_rank < i_rank) {
                            int rc = REDUCED_COST( i, j, a );
                            int j_new_rank;
                            if (rc < 0) {
                                j_new_rank = i_rank;
                            }else{
                                int dr = rc / epsilon;
                                j_new_rank = ( dr < linf ) ? i_rank - ( dr + 1 ) : 0;
                            }//if
                            if (j_rank < j_new_rank) {
                                if (cc == 1) {
                                    j->rank = j_new_rank;
                                    if ( j_rank > 0 ) {
                                        McmfBucket *b_old = buckets + j_rank;
                                        remove_from_bucket( j, b_old );
                                    }
                                    McmfBucket *b_new = buckets + j_new_rank;
                                    insert_to_bucket( j, b_new );
                                }else {
                                    int df = a->rez_capacity;
                                    increase_flow( i, j, a, df );
                                }//if
                            }//if
                        }//if
                    }//if
                    a++;
                }//while
                i->price -= dp;
            }//while
        }//for
        if (cc == 0)
            break;
    }//while
    
    if (cc == 0) {
        for (McmfNode * i = nodes; i != sentinel_node; i ++) {
            McmfEdge *a = i->first;
            McmfEdge *a_stop = (i + 1)->suspended;
            while (a != a_stop) {
                if ( REDUCED_COST( i, a->head, a ) < - epsilon ) {
                    int df = a->rez_capacity;
                    if (df > 0) {
                        increase_flow( i, a->head, a, df );
                    }
                }
                a++;
            }
        }
    }
    
    return cc;
}

void Mcmf::price_out(){
    double n_cut_off = - cut_off;
    for ( McmfNode *i = nodes; i != sentinel_node; i ++) {
        McmfEdge *a = i->first;
        McmfEdge *a_stop = (i + 1)->suspended;
        while (a != a_stop) {
            double rc = REDUCED_COST( i, a->head, a );
            if ( ( rc > cut_off && a->sister->rez_capacity <= 0 ) ||
                ( rc < n_cut_off && a->rez_capacity <= 0 ) ) { // 终止边
                
                McmfEdge *b = i->first;
                i->first++;
                exchange(a, b);
            }
            a++;
        }
    }
}


// 建立结果
void Mcmf::build_result( double *objective_cost){
    McmfEdge *a = arcs;
    int na = 0;
    double obj_internal = 0;
    
    while(a != sentinel_arc)
    {
        int flow = cap[na] - a->rez_capacity;
        if(cap[na] > 0 && flow != 0 )
        {
            int cs = a->cost / dn;
            obj_internal += (double) cs * (double) flow;
            //printf("%d-%d:%d %d\n", a->_tail- _nodes, a->_head - _nodes, flow, cs);
            int tail_id = (int)(a->tail - nodes);
            int head_id = (int)(a->head - nodes);
            this->roads.push_back(Edge(tail_id, head_id, flow, cs));
            if (tail_id == n-2){
                this->server_supply[head_id] = flow;
            }
        }
        
        a++;
        na++;
    }
    
    *objective_cost = obj_internal;
}

// 核心算法
void Mcmf::algo(double *objective_cost){
    int cc = 0;
    update_epsilon();
    
    do {
        refine();
        if (this -> is_error){
            return;
        }
        if (n_ref >= 1){
            price_out();
        }
        if (update_epsilon()){
            break;
        }
        while (1) {
            if (!price_refine()) {
                break;
            }
            if ( n_ref >= 1 ) {
                if (price_in()){
                    break;
                }
                if ((cc = update_epsilon())){
                    break;
                }
            }
        }
    } while ( cc == 0 );
    
    if (this->is_error){
        return;
    }
    
    build_result(objective_cost);
}

// 执行
int Mcmf::exec() {
    
    double objective_cost = 0;
    
    pre_processing();
    
    m = 2 * m;
    algo_initialize();
    
    algo(&objective_cost);
    
    //	PRINT("cost:%.2f\n", objective_cost);
    
    if(!this -> is_error)
    {
        this->cost_maxflow = objective_cost ;
    }
   
    //	printf("McmfEdge:%ld, fail:%ld\n", this->roads.size(), this->demand_situation.size() );
    
    deallocate_arrays();
    return 0;
}



// 共享utils
void Mcmf::increase_flow( McmfNode *i, McmfNode *j, McmfEdge *a, long df) {
    i->excess -= df;
    j->excess += df;
    a->rez_capacity -= df;
    a->sister->rez_capacity += df;
}


// 超量队列的utils
void Mcmf::reset_excess_q() {
    for ( ; excq_first != NULL; excq_first = excq_last ) {
        excq_last = excq_first->q_next;
        excq_first->q_next = sentinel_node;
    }
}

void Mcmf::insert_to_excess_q( McmfNode *i) {
    if ( excq_first != NULL ) {
        excq_last->q_next = i;
    } else {
        excq_first = i;
    }
    i->q_next = NULL;
    excq_last = i;
}

void Mcmf::insert_to_bucket( McmfNode *i, McmfBucket *b) {
    i->b_next = b->p_first;
    b->p_first->b_prev = i;
    b->p_first = i;
}

void Mcmf::remove_from_bucket( McmfNode *i, McmfBucket *b) {
    if ( i == b->p_first ) {
        b->p_first = i->b_next;
    } else {
        i->b_prev->b_next = i->b_next;
        i->b_next->b_prev = i->b_prev;
    }
}

void Mcmf::update_cut_off() {
    if ( n_bad_pricein + n_bad_relabel == 0) {
        cut_off_factor = 1 * pow( (double)n, 0.75 );
        cut_off_factor = MAX ( cut_off_factor, 12 );
        cut_off = cut_off_factor * epsilon;
        cut_on = cut_off * 0.8;
    } else {
        cut_off_factor *= 4;
        cut_off = cut_off_factor * epsilon;
        cut_on = cut_off * 0.8;
    }
}
void Mcmf::exchange( McmfEdge *a, McmfEdge *b) {
    if ( a != b) {
        McmfEdge *sa = a->sister;
        McmfEdge *sb = b->sister;
        int d_cap;
        
        d_arc.rez_capacity = a->rez_capacity;
        d_arc.cost = a->cost;
        d_arc.head = a->head;
        d_arc.tail = a->tail;
        
        a->rez_capacity = b->rez_capacity;
        a->cost = b->cost;
        a->head = b->head;
        a->tail = b->tail;
        
        b->rez_capacity = d_arc.rez_capacity;
        b->cost = d_arc.cost;
        b->head = d_arc.head;
        b->tail = d_arc.tail;
        
        if ( a != sb) {
            b->sister = sa;
            a->sister = sb;
            sa->sister = b;
            sb->sister = a;
        }
        
        d_cap = cap[ a - arcs];
        cap[ a - arcs] = cap[ b - arcs];
        cap[ b - arcs] = d_cap;
    }
}

void MfGraph::initGraph()
{
    initGraphFast();
    topLevelS = topLevelT = 1;
}


void MfGraph::initSize(int numNodes, int numEdges)
{
    unsigned long long arcTmpMemsize = (unsigned long long)sizeof(TmpEdge)*(unsigned long long)numEdges;
    unsigned long long arcRealMemsize = (unsigned long long)sizeof(Arc)*(unsigned long long)(numEdges*2);
    unsigned long long nodeMemsize = (unsigned long long)sizeof(Node**)*(unsigned long long)(numNodes*3) +
    (IB_EXCESSES ? ((unsigned long long)sizeof(Node**)*(unsigned long long)(numNodes*2)) : 0);
    unsigned long long arcMemsize = 0;

    arcMemsize = arcRealMemsize + arcTmpMemsize;

    if (arcMemsize < (arcRealMemsize + nodeMemsize)) {
        arcMemsize = (arcRealMemsize + nodeMemsize);
    }

    memArcs = new char[arcMemsize];
    memset(memArcs, 0, (unsigned long long)sizeof(char)*arcMemsize);

    tmpEdges = (TmpEdge*)(memArcs + arcRealMemsize);

    tmpEdgeLast = tmpEdges;
    arcs = (Arc*)memArcs;
    arcEnd = arcs + numEdges*2;

    this->numNodes = numNodes;
    nodes = new Node[numNodes+1];
    memset(nodes, 0, sizeof(Node)*(numNodes+1));
    nodeEnd = nodes+numNodes;
    active0.init((Node**)(arcEnd));
    activeS1.init((Node**)(arcEnd) + numNodes);
    activeT1.init((Node**)(arcEnd) + (2*numNodes));
    if (IB_EXCESSES) {
        ptrs = (Node**)(arcEnd) + (3*numNodes);
        excessBuckets.init(nodes, ptrs, numNodes);
    }
    orphan3PassBuckets.init(nodes, numNodes);
    orphanBuckets.init(nodes, numNodes);

    flow = 0;
}


void MfGraph::initNodes()
{
    Node *x;
    for (x=nodes; x <= nodeEnd; x++) {
        x->firstArc = (arcs + x->label);
        if (x->excess == 0) {
            x->label = 0;
            continue;
        }
        if (x->excess > 0) {
            x->label = 1;
            activeS1.add(x);
        } else {
            x->label = -1;
            activeT1.add(x);
        }
    }
}

void MfGraph::initGraphFast()
{
    Node *x;
    TmpEdge *te;
    Arc *a;

    nodes->firstArc = arcs;
    for (x=nodes; x != nodeEnd; x++) {
        (x+1)->firstArc = x->firstArc + x->label;
        x->label = (int)(x->firstArc-arcs);
    }
    nodeEnd->label = (int)(arcEnd-arcs);

    for (te=tmpEdges; te != tmpEdgeLast; te++) {
        a = (nodes+te->tail)->firstArc;
        a->rev = (nodes+te->head)->firstArc;
        a->head = nodes+te->head;
        a->rCap = te->cap;
        a->isRevResidual = (te->revCap != 0);

        a = (nodes+te->head)->firstArc;
        a->rev = (nodes+te->tail)->firstArc;
        a->head = nodes+te->tail;
        a->rCap = te->revCap;
        a->isRevResidual = (te->cap != 0);

        ++((nodes+te->head)->firstArc);
        ++((nodes+te->tail)->firstArc);
    }

    initNodes();
}

template<bool sTree> int MfGraph::augmentPath(Node *x, int push)
{
    Node *y;
    Arc *a;
    int orphanMinLevel = (sTree ? topLevelS : topLevelT) + 1;

    augTimestamp++;
    for (; ; x=a->head)
    {
        if (x->excess) break;
        a = x->parent;
        if (sTree) {
            a->rCap += push;
            a->rev->isRevResidual = 1;
            a->rev->rCap -= push;
        } else {
            a->rev->rCap += push;
            a->isRevResidual = 1;
            a->rCap -= push;
        }

        if ((sTree ? (a->rev->rCap) : (a->rCap)) == 0)
        {
            if (sTree) a->isRevResidual = 0;
            else a->rev->isRevResidual = 0;
            REMOVE_SIBLING(x,y);
            orphanMinLevel = (sTree ? x->label : -x->label);
            orphanBuckets.add<sTree>(x);
        }
    }
    x->excess += (sTree ? -push : push);
    if (x->excess == 0) {
        orphanMinLevel = (sTree ? x->label : -x->label);
        orphanBuckets.add<sTree>(x);
    }
    flow += push;

    return orphanMinLevel;
}


template<bool sTree> int MfGraph::augmentExcess(Node *x, int push)
{
    Node *y;
    Arc *a;
    int orphanMinLevel = (sTree ? topLevelS : topLevelT)+1;
    augTimestamp++;
    while (sTree ? (x->excess <= 0) : (x->excess >= 0))
    {

        a = x->parent;

        if ((sTree ? (a->rev->rCap) : (a->rCap)) < (sTree ? (push-x->excess) : (x->excess+push))) {
            x->excess += (sTree ? (a->rev->rCap - push) : (push-a->rCap));
            push = (sTree ? a->rev->rCap : a->rCap);
        } else {
            push += (sTree ? -(x->excess) : x->excess);
            x->excess = 0;
        }

        if (sTree) {
            a->rCap += push;
            a->rev->isRevResidual = 1;
            a->rev->rCap -= push;
        } else {
            a->rev->rCap += push;
            a->isRevResidual = 1;
            a->rCap -= push;
        }

        if ((sTree ? (a->rev->rCap) : (a->rCap)) == 0)
        {
            if (sTree) a->isRevResidual = 0;
            else a->rev->isRevResidual = 0;
            REMOVE_SIBLING(x,y);
            orphanMinLevel = (sTree ? x->label : -x->label);
            orphanBuckets.add<sTree>(x);
            if (x->excess) excessBuckets.incMaxBucket(sTree ? x->label : -x->label);
        }

        x = a->head;
        if (sTree ? (x->excess < 0) : (x->excess > 0)) excessBuckets.remove<sTree>(x);
    }

    if (push <= (sTree ? (x->excess) : -(x->excess))) flow += push;
    else flow += (sTree ? (x->excess) : -(x->excess));
    x->excess += (sTree ? (-push) : push);
    if (sTree ? (x->excess <= 0) : (x->excess >= 0)) {
        orphanMinLevel = (sTree ? x->label : -x->label);
        orphanBuckets.add<sTree>(x);
        if (x->excess) excessBuckets.incMaxBucket(sTree ? x->label : -x->label);
    }

    return orphanMinLevel;
}

template<bool sTree> void MfGraph::augmentExcesses()
{
    Node *x;
    int minOrphanLevel;
    int adoptedUpToLevel = excessBuckets.maxBucket;

    if (!excessBuckets.empty())
        for (; excessBuckets.maxBucket != (excessBuckets.minBucket-1); excessBuckets.maxBucket--)
            while ((x=excessBuckets.popFront(excessBuckets.maxBucket)) != NULL)
            {
                minOrphanLevel = augmentExcess<sTree>(x, 0);
                if (adoptedUpToLevel < minOrphanLevel) minOrphanLevel = adoptedUpToLevel;
                adoption<sTree>(minOrphanLevel, false);
                adoptedUpToLevel = excessBuckets.maxBucket;
            }
    excessBuckets.reset();
    if (orphanBuckets.maxBucket != 0) adoption<sTree>(adoptedUpToLevel+1, true);
    while ((x=excessBuckets.popFront(0)) != NULL) orphanFree<sTree>(x);
}

void MfGraph::augment(Arc *bridge)
{
    Node *x, *y;
    Arc *a;
    int bottleneck, bottleneckT, bottleneckS, minOrphanLevel;
    bool forceBottleneck;

    forceBottleneck = (IB_EXCESSES ? false : true);
    if (IB_BOTTLENECK_ORIG && IB_EXCESSES)
    {
        bottleneck = bridge->rCap;
        if (bridge->head->excess != 0 && -(bridge->head->excess) < bottleneck) {
            bottleneck = -(bridge->head->excess);
        }
        if (bridge->rev->head->excess != 0 && bridge->rev->head->excess < bottleneck) {
            bottleneck = bridge->rev->head->excess;
        }
    }
    else
    {
        bottleneck = bottleneckS = bridge->rCap;
        if (bottleneck != 1) {
            for (x=bridge->rev->head; ; x=a->head)
            {
                if (x->excess) break;
                a = x->parent;
                if (bottleneckS > a->rev->rCap) {
                    bottleneckS = a->rev->rCap;
                }
            }
            if (bottleneckS > x->excess) {
                bottleneckS = x->excess;
            }
            if (IB_EXCESSES && x->label != 1) forceBottleneck = true;
            if (x == bridge->rev->head) bottleneck = bottleneckS;
        }

        if (bottleneck != 1) {
            bottleneckT = bridge->rCap;
            for (x=bridge->head; ; x=a->head)
            {
                if (x->excess) break;
                a = x->parent;
                if (bottleneckT > a->rCap) {
                    bottleneckT = a->rCap;
                }
            }
            if (bottleneckT > (-x->excess)) {
                bottleneckT = (-x->excess);
            }
            if (IB_EXCESSES && x->label != -1) forceBottleneck = true;
            if (x == bridge->head && bottleneck > bottleneckT) bottleneck = bottleneckT;

            if (forceBottleneck) {
                if (bottleneckS < bottleneckT) bottleneck = bottleneckS;
                else bottleneck = bottleneckT;
            }
        }
    }

    if (IBSTATS) {
        int augLen = (-(bridge->head->label)-1 + bridge->rev->head->label-1 + 1);
    }

    bridge->rev->rCap += bottleneck;
    bridge->isRevResidual = 1;
    bridge->rCap -= bottleneck;
    if (bridge->rCap == 0) {
        bridge->rev->isRevResidual = 0;
    }
    flow -= bottleneck;

    x = bridge->head;
    if (!IB_EXCESSES || bottleneck == 1 || forceBottleneck) {
        minOrphanLevel = augmentPath<false>(x, bottleneck);
        adoption<false>(minOrphanLevel, true);
    } else if (IB_ADOPTION_PR && !x->excess) {
        x->excess += bottleneck;
        excessBuckets.add<false>(x);
        REMOVE_SIBLING(x,y);
        augmentExcessesDischarge<false>();
    } else {
        minOrphanLevel = augmentExcess<false>(x, bottleneck);
        adoption<false>(minOrphanLevel, false);
        augmentExcesses<false>();
    }

    x = bridge->rev->head;
    if (!IB_EXCESSES || bottleneck == 1 || forceBottleneck) {
        minOrphanLevel = augmentPath<true>(x, bottleneck);
        adoption<true>(minOrphanLevel, true);
    } else if (IB_ADOPTION_PR && !x->excess) {
        x->excess -= bottleneck;
        excessBuckets.add<true>(x);
        REMOVE_SIBLING(x,y);
        augmentExcessesDischarge<true>();
    } else {
        minOrphanLevel = augmentExcess<true>(x, bottleneck);
        adoption<true>(minOrphanLevel, false);
        augmentExcesses<true>();
    }
}

template<bool sTree> void MfGraph::adoption(int fromLevel, bool toTop)
{
    Node *x, *y, *z;
    Arc *a;
    Arc *aEnd;
    int threePassLevel;
    int minLabel, numOrphans, numOrphansUniq;
    int level;

    threePassLevel=0;
    numOrphans=0;
    numOrphansUniq=0;
    for (level = fromLevel;
         level <= orphanBuckets.maxBucket &&
         (!IB_EXCESSES || toTop || threePassLevel || level <= excessBuckets.maxBucket);
         level++)
        while ((x=orphanBuckets.popFront(level)) != NULL)
        {
            numOrphans++;
            if (x->lastAugTimestamp != augTimestamp) {
                x->lastAugTimestamp = augTimestamp;
                if (sTree) uniqOrphansS++;
                else uniqOrphansT++;
                numOrphansUniq++;
            }
            if (IB_HYBRID_ADOPTION && threePassLevel == 0 && numOrphans >= 3*numOrphansUniq) {
                threePassLevel = 1;
            }

            if (x->isParentCurr) {
                a = x->parent;
            } else {
                a = x->firstArc;
                x->isParentCurr = 1;
            }
            x->parent = NULL;
            aEnd = (x+1)->firstArc;
            if (x->label != (sTree ? 1 : -1))
            {
                minLabel = x->label - (sTree ? 1 : -1);
                for (; a != aEnd; a++)
                {
                    y = a->head;
                    if ((sTree ? a->isRevResidual : a->rCap) != 0 && y->label == minLabel)
                    {
                        x->parent = a;
                        ADD_SIBLING(x,y);
                        break;
                    }
                }
            }
            if (x->parent != NULL) {
                if (IB_EXCESSES && x->excess) excessBuckets.add<sTree>(x);
                continue;
            }


            if (x->label == (sTree ? topLevelS : -topLevelT)) {
                orphanFree<sTree>(x);
                continue;
            }
            for (y=x->firstSon; y != NULL; y=z)
            {
                z=y->nextPtr;
                if (IB_EXCESSES && y->excess) excessBuckets.remove<sTree>(y);
                orphanBuckets.add<sTree>(y);
            }
            x->firstSon = NULL;

            if (threePassLevel) {
                x->label += (sTree ? 1 : -1);
                orphan3PassBuckets.add<sTree>(x);
                if (threePassLevel == 1) {
                    threePassLevel = level+1;
                }
                continue;
            }

            minLabel = (sTree ? topLevelS : -topLevelT);
            if (x->label != minLabel) for (a=x->firstArc; a != aEnd; a++)
            {
                y = a->head;
                if ((sTree ? a->isRevResidual : a->rCap) &&
                    (sTree ? (y->label > 0) : (y->label < 0)) &&
                    (sTree ? (y->label < minLabel) : (y->label > minLabel)))
                {
                    minLabel = y->label;
                    x->parent = a;
                    if (minLabel == x->label) break;
                }
            }

            if (x->parent != NULL) {
                x->label = minLabel + (sTree ? 1 : -1);
                ADD_SIBLING(x, x->parent->head);
                if (sTree) {
                    if (x->label == topLevelS) activeS1.add(x);
                } else {
                    if (x->label == -topLevelT) activeT1.add(x);
                }
                if (IB_EXCESSES && x->excess) excessBuckets.add<sTree>(x);
            } else {
                orphanFree<sTree>(x);
            }
        }
    if (level > orphanBuckets.maxBucket) orphanBuckets.maxBucket=0;

    if (threePassLevel) {
        adoption3Pass<sTree>(threePassLevel);
    }
}

template <bool sTree> void MfGraph::adoption3Pass(int minBucket)
{
    Arc *a, *aEnd;
    Node *x, *y;
    int minLabel, destLabel;

    for (int level=minBucket; level <= orphan3PassBuckets.maxBucket; level++)
        while ((x = orphan3PassBuckets.popFront(level)) != NULL)
        {
            aEnd = (x+1)->firstArc;

            if (x->parent == NULL) {
                minLabel = (sTree ? topLevelS : -topLevelT);
                destLabel = x->label - (sTree ? 1 : -1);
                for (a=x->firstArc; a != aEnd; a++) {
                    y = a->head;
                    if ((sTree ? a->isRevResidual : a->rCap) &&
                        ((sTree ? (y->excess > 0) : (y->excess < 0)) || y->parent != NULL) &&
                        (sTree ? (y->label > 0) : (y->label < 0)) &&
                        (sTree ? (y->label < minLabel) : (y->label > minLabel)))
                    {
                        x->parent = a;
                        if ((minLabel = y->label) == destLabel) break;
                    }
                }
                if (x->parent == NULL) {
                    x->label = 0;
                    if (IB_EXCESSES && x->excess) excessBuckets.add<sTree>(x);
                    continue;
                }
                x->label = minLabel + (sTree ? 1 : -1);
                if (x->label != (sTree ? level : -level)) {
                    orphan3PassBuckets.add<sTree>(x);
                    continue;
                }
            }

            if (x->label != (sTree ? topLevelS : -topLevelT))
            {
                minLabel = x->label + (sTree ? 1 : -1);
                for (a=x->firstArc; a != aEnd; a++) {
                    y = a->head;

                    if ((sTree ? a->rCap : a->isRevResidual) &&
                        (y->label == 0 ||
                         (sTree ? (minLabel < y->label) : (minLabel > y->label))))
                    {
                        if (y->label != 0) orphan3PassBuckets.remove<sTree>(y);
                        else if (IB_EXCESSES && y->excess) excessBuckets.remove<sTree>(y);
                        y->label = minLabel;
                        y->parent = a->rev;
                        orphan3PassBuckets.add<sTree>(y);
                    }
                }
            }

            ADD_SIBLING(x, x->parent->head);
            x->isParentCurr = 0;
            if (IB_EXCESSES && x->excess) excessBuckets.add<sTree>(x);

            if (sTree) {
                if (x->label == topLevelS) activeS1.add(x);
            } else {
                if (x->label == -topLevelT) activeT1.add(x);
            }
        }
    orphan3PassBuckets.maxBucket = 0;
}


template<bool dirS> void MfGraph::growth()
{
    Node *x, *y;
    Arc *a, *aEnd;

    for (Node **active=active0.list; active != (active0.list + active0.len); active++)
    {
        x = (*active);

        if (x->label != (dirS ? (topLevelS-1): -(topLevelT-1))) {
            continue;
        }

        aEnd = (x+1)->firstArc;
        for (a=x->firstArc; a != aEnd; a++)
        {
            if ((dirS ? a->rCap : a->isRevResidual) == 0) continue;
            y = a->head;
            if (y->label == 0)
            {
                y->isParentCurr = 0;
                y->label = x->label + (dirS ? 1 : -1);
                y->parent = a->rev;
                ADD_SIBLING(y, x);
                if (dirS) activeS1.add(y);
                else activeT1.add(y);
            }
            else if (dirS ? (y->label < 0) : (y->label > 0))
            {
                augment(dirS ? a : (a->rev));
                if (x->label != (dirS ? (topLevelS-1) : -(topLevelT-1))) {
                    break;
                }
                if (dirS ? (a->rCap) : (a->isRevResidual)) a--;
            }
        }
    }
    active0.clear();
}

template<bool sTree> void MfGraph::augmentIncrements()
{
    Node *x, *y;
    Node **end = incList+incLen;
    int minOrphanLevel = 1<<30;

    for (Node **inc=incList; inc != end; inc++)
    {
        x = (*inc);
        if (!x->isIncremental || (sTree ? (x->label < 0) : (x->label > 0))) continue;
        x->isIncremental = 0;
        if (x->label == 0)
        {
            if (!x->excess) continue;
            x->isParentCurr = 0;
            if (x->excess > 0) {
                x->label = topLevelS;
                activeS1.add(x);
            } else if (x->excess < 0) {
                x->label = -topLevelT;
                activeT1.add(x);
            }
        }
        else if ((sTree ? (x->excess <= 0) : (x->excess >= 0)) &&
                 (!x->parent || !(sTree ? x->parent->isRevResidual : x->parent->rCap)))
        {
            if (x->parent) REMOVE_SIBLING(x,y);
            if ((sTree ? x->label : -x->label) < minOrphanLevel) {
                minOrphanLevel = (sTree ? x->label : -x->label);
            }
            orphanBuckets.add<sTree>(x);
            if (x->excess) excessBuckets.incMaxBucket(sTree ? x->label : -x->label);
        }
        else if (sTree ? (x->excess < 0) : (x->excess > 0))
        {
            excessBuckets.add<sTree>(x);
        }
        else if (x->excess && x->parent)
        {
            REMOVE_SIBLING(x,y);
            x->parent = NULL;
            x->isParentCurr = 0;
        }
    }
    if (orphanBuckets.maxBucket != 0) adoption<sTree>(minOrphanLevel, false);
    if (IB_ADOPTION_PR) augmentExcessesDischarge<sTree>();
    else augmentExcesses<sTree>();
}

int MfGraph::computeMaxFlow()
{
    return computeMaxFlow(true, false);
}

int MfGraph::computeMaxFlow(bool allowIncrements)
{
    return computeMaxFlow(true, allowIncrements);
}

int MfGraph::computeMaxFlow(bool initialDirS, bool allowIncrements)
{
    if (incIteration >= 1 && incList != NULL) {
        augmentIncrements<true>();
        augmentIncrements<false>();
        incList = NULL;
    }
    bool dirS = initialDirS;
    while (true)
    {
        if (dirS) {
            ActiveList::swapLists(&active0, &activeS1);
            topLevelS++;
        } else {
            ActiveList::swapLists(&active0, &activeT1);
            topLevelT++;
        }
        orphanBuckets.allocate((topLevelS > topLevelT) ? topLevelS : topLevelT);
        orphan3PassBuckets.allocate((topLevelS > topLevelT) ? topLevelS : topLevelT);
        if (IB_EXCESSES) excessBuckets.allocate((topLevelS > topLevelT) ? topLevelS : topLevelT);
        if (dirS) growth<true>();
        else growth<false>();

        if (!allowIncrements && (activeS1.len == 0 || activeT1.len == 0)) break;
        if (activeS1.len == 0 && activeT1.len == 0) break;
        if (activeT1.len == 0) dirS=true;
        else if (activeS1.len == 0) dirS=false;
        else if (!IB_ALTERNATE_SMART && dirS) dirS = false;
        else if (IB_ALTERNATE_SMART && uniqOrphansT == uniqOrphansS && dirS) dirS=false;
        else if (IB_ALTERNATE_SMART && uniqOrphansT < uniqOrphansS) dirS=false;
        else dirS=true;
    }

    incIteration++;
    return flow;
}

template<bool sTree> void MfGraph::augmentExcessesDischarge()
{
    Node *x;
    if (!excessBuckets.empty())
        for (; excessBuckets.maxBucket != (excessBuckets.minBucket-1); excessBuckets.maxBucket--)
            while ((x=excessBuckets.popFront(excessBuckets.maxBucket)) != NULL) {
                augmentDischarge<sTree>(x);
            }
    excessBuckets.reset();
    while ((x=excessBuckets.popFront(0)) != NULL) {
        x->isIncremental = 0;
        orphanBuckets.add<sTree>(x);
    }
    augTimestamp++;
    adoption<sTree>(1, true);
}

template<bool sTree> void MfGraph::augmentDischarge(Node *x)
{
    Node *y, *z;
    int minLabel, push;
    Arc *aEnd = (x+1)->firstArc;
    Arc *a;
    int startLabel = x->label;

    while (true)
    {
        if (x->isParentCurr) {
            a = x->parent;
        } else {
            a = x->firstArc;
            x->isParentCurr = 1;
        }
        if (x->label != (sTree ? 1 : -1) && a != NULL)
        {
            minLabel = x->label + (sTree ? -1 : 1);
            for (; a != aEnd; a++)
            {
                y = a->head;
                if ((sTree ? a->isRevResidual : a->rCap) == 0 || y->label != minLabel) {
                    continue;
                }

                push = (sTree ? (a->rev->rCap) : (a->rCap));
                if (push > (sTree ? (-x->excess) : (x->excess))) {
                    push = (sTree ? (-x->excess) : (x->excess));
                }
                x->excess += (sTree ? push : (-push));
                if (sTree) {
                    a->rev->rCap -= push;
                    a->rCap += push;
                    a->rev->isRevResidual = 1;
                    a->isRevResidual = (a->rev->rCap ? 1 : 0);
                } else {
                    a->rCap -= push;
                    a->rev->rCap += push;
                    a->rev->isRevResidual = (a->rCap ? 1 : 0);
                    a->isRevResidual = 1;
                }

                if (sTree && y->excess > 0) {
                    if (y->excess >= push) flow += push;
                    else flow += y->excess;
                } else if (!sTree && y->excess < 0) {
                    if (-y->excess >= push) flow += push;
                    else flow -= y->excess;
                }
                y->excess += (sTree ? (-push) : push);
                if (y->excess == 0 ) {
                    y->label = 0;
                    excessBuckets.add<sTree>(y);
                    y->label = minLabel;
                    y->isIncremental = 1;
                } else if (sTree ? (y->excess < 0 && y->excess >= -push) : (y->excess > 0 && y->excess <= push)) {
                    if (y->isIncremental) {
                        y->label = 0;
                        excessBuckets.remove<sTree>(y);
                        y->label = minLabel;
                        y->isIncremental = 0;
                    } else if (y->parent != NULL) {
                        REMOVE_SIBLING(y,z);
                    }
                    excessBuckets.add<sTree>(y);
                }
                if (x->excess == 0) {
                    x->parent = a;
                    if (!(sTree ? a->isRevResidual : a->rCap)) {
                        x->label = 0;
                        excessBuckets.add<sTree>(x);
                        x->label = minLabel + (sTree ? 1 : -1);
                        x->isIncremental = 1;
                    }
                    break;
                }
            }
        }
        if (x->excess == 0) break;

        minLabel = x->label + (sTree ? 1 : -1);
        for (y=x->firstSon; y != NULL; y=z)
        {
            z=y->nextPtr;
            y->label = 0;
            excessBuckets.add<sTree>(y);
            y->label = minLabel;
            y->isIncremental = 1;
        }
        x->firstSon = NULL;
        
        minLabel = (sTree ? topLevelS : -topLevelT);
        x->parent = NULL;
        for (a=x->firstArc; a != aEnd; a++)
        {
            y = a->head;
            if ((sTree ? a->isRevResidual : a->rCap) &&
                (sTree ? (y->label > 0) : (y->label < 0)) &&
                (sTree ? (y->label < minLabel) : (y->label > minLabel)))
            {
                minLabel = y->label;
                x->parent = a;
                if (minLabel == x->label) break;
            }
        }
        if (x->parent != NULL) {
            x->label = minLabel + (sTree ? 1 : -1);
        } else {
            orphanFree<sTree>(x);
            break;
        }
    }
    
    if (x->parent != NULL && !x->isIncremental) ADD_SIBLING(x, x->parent->head);
    if (sTree) {
        if (startLabel != x->label && x->label == topLevelS) activeS1.add(x);
    } else {
        if (startLabel != x->label && x->label == -topLevelT) activeT1.add(x);
    }
}

void MfGraph::pushRelabelShelve(int fromLevel)
{
    Node *x = NULL;
    for (int bucket=fromLevel; bucket <= prNodeBuckets.maxBucket; bucket++) {
        if (prNodeBuckets.isEmpty(bucket)) continue;
        
        for (x=prNodeBuckets.buckets[bucket]; x != NULL; x = x->nextPtr) x->label = 0;
    }
    int numLevels = prNodeBuckets.maxBucket - fromLevel + 1;
    memset(prNodeBuckets.buckets + fromLevel, 0, sizeof(Node*)*numLevels);
    memset(excessBuckets.buckets + fromLevel, 0, sizeof(Node*)*numLevels);
    prNodeBuckets.maxBucket = fromLevel-1;
    excessBuckets.maxBucket = fromLevel-1;
}

void MfGraph::pushRelabel()
{
    return pushRelabelDir<false>();
}

template<bool sTree> void MfGraph::pushRelabelDir()
{
    Node *x;
    int level;
    
    topLevelS = topLevelT = numNodes;
    pushRelabelGlobalUpdate<sTree>();
    
    int nDischarges = 0;
    for (; excessBuckets.maxBucket >= excessBuckets.minBucket; excessBuckets.maxBucket--)
        while ((x=excessBuckets.popFront(excessBuckets.maxBucket)) != NULL)
        {
            
            level = excessBuckets.maxBucket;
            pushRelabelDischarge<sTree>(x);
            nDischarges++;
            if (prNodeBuckets.maxBucket < level) {
                excessBuckets.allocate(level+2);
                prNodeBuckets.allocate(level+2);
            }
            
            if (nDischarges % (30*numNodes) == 0) pushRelabelGlobalUpdate<sTree>();
            else if (prNodeBuckets.isEmpty(level)) pushRelabelShelve(level+1);
        }
}

template<bool sTree> void MfGraph::pushRelabelGlobalUpdate()
{
    Node *x, *y;
    Arc *a, *aEnd;
    
    memset(prNodeBuckets.buckets, 0, sizeof(Node*)*(prNodeBuckets.allocLevels+1));
    memset(excessBuckets.buckets, 0, sizeof(Node*)*(excessBuckets.allocLevels+1));
    prNodeBuckets.maxBucket = 1;
    excessBuckets.reset();
    for (x=nodes; x != nodeEnd; x++) {
        x->parent = NULL;
        x->isParentCurr = 0;
        if ((sTree ? (x->excess > 0) : (x->excess < 0))) {
            x->label = (sTree ? 1 : -1);
            prNodeBuckets.add<sTree>(x);
        } else x->label = 0;
    }
    for (int bucket=1; bucket <= prNodeBuckets.maxBucket; bucket++)
        for (x=prNodeBuckets.buckets[bucket]; x != NULL; x = x->nextPtr) {
            aEnd = (x+1)->firstArc;
            for (a=x->firstArc; a != aEnd; a++) {
                if (!(sTree ? a->rCap : a->isRevResidual)) continue;
                y = a->head;
                if (y->parent != NULL || (sTree ? (y->excess > 0) : (y->excess < 0))) continue;
                y->label = (sTree ? (bucket+1) : (-bucket-1));
                prNodeBuckets.add<sTree>(y);
                y->parent = a->rev;
                if (y->excess) excessBuckets.add<sTree>(y);
            }
        }
}

template<bool sTree> void MfGraph::pushRelabelDischarge(Node *x)
{
    Node *y;
    int minLabel, push;
    Arc *aEnd = (x+1)->firstArc;
    Arc *a;
    
    prNodeBuckets.remove<sTree>(x);
    while (true)
    {
        if (x->isParentCurr) {
            a = x->parent;
        } else {
            a = x->firstArc;
            x->isParentCurr = 1;
        }
        if (x->label != (sTree ? 1 : -1))
        {
            minLabel = x->label - (sTree ? 1 : -1);
            for (; a != aEnd; a++)
            {
                y = a->head;
                if ((sTree ? a->isRevResidual : a->rCap) == 0 || y->label != minLabel) {
                    continue;
                }
                
                push = (sTree ? (a->rev->rCap) : (a->rCap));
                if (push > (sTree ? (-x->excess) : (x->excess))) {
                    push = (sTree ? (-x->excess) : (x->excess));
                }
                x->excess += (sTree ? push : (-push));
                if (sTree) {
                    a->rev->rCap -= push;
                    a->rCap += push;
                    a->rev->isRevResidual = 1;
                    a->isRevResidual = (a->rev->rCap ? 1 : 0);
                } else {
                    a->rCap -= push;
                    a->rev->rCap += push;
                    a->rev->isRevResidual = (a->rCap ? 1 : 0);
                    a->isRevResidual = 1;
                }
                
                if (sTree && y->excess > 0) {
                    if (y->excess >= push) flow += push;
                    else flow += y->excess;
                } else if (!sTree && y->excess < 0) {
                    if (-y->excess >= push) flow += push;
                    else flow -= y->excess;
                }
                y->excess += (sTree ? (-push) : push);
                if (sTree ? (y->excess < 0 && y->excess >= -push) : (y->excess > 0 && y->excess <= push)) {
                    excessBuckets.add<sTree>(y);
                }
                if (x->excess == 0) {
                    x->parent = a;
                    break;
                }
            }
        }
        if (x->excess == 0) break;
        
        minLabel = (sTree ? (numNodes-1) : (-numNodes+1));
        x->parent = NULL;
        for (a=x->firstArc; a != aEnd; a++)
        {
            y = a->head;
            if ((sTree ? a->isRevResidual : a->rCap) &&
                (sTree ? (y->label > 0) : (y->label < 0)) &&
                (sTree ? (y->label < minLabel) : (y->label > minLabel)))
            {
                minLabel = y->label;
                x->parent = a;
                if (minLabel == x->label) break;
            }
        }
        if (x->parent != NULL) {
            x->label = minLabel + (sTree ? 1 : -1);
        } else {
            x->label = 0;
            break;
        }
    }
    if (x->label != 0) prNodeBuckets.add<sTree>(x);
}

void MfGraph::analysis_fail(){
    Node *sink = nodes + numNodes - 1;
    Arc *a;
    
    for (a = sink->firstArc; a != (sink+1)->firstArc; a++){
        if (a->rev->rCap != 0) {
            this->fail_consumer.push_back(a->head - nodes);
        }
    }
}
