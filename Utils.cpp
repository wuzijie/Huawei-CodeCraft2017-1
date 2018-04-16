//
//  Utils.cpp
//  CDN复赛
//
//  Created by Daniel on 2017/4/9.
//  Copyright © 2017年 Daniel. All rights reserved.
//


#include "cdn.h"

/**
 时间检查
 
 @return 超时返回ture
 */
bool TimeOut(){
    struct timeb time;
    ftime(&time);
    
    static unsigned long start1 = time.time*1000 + time.millitm;
    unsigned long end = time.time*1000 + time.millitm;
    
    if (end - start1 > TIME_OUT*1000) {
        return true;
    }
    return false;
}
