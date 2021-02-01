//
//  structs.h
//  ShaPMoD
//
//  Created by Alessandro Coretti on 12/15/17.
//  Copyright © 2017 Alessandro Coretti. All rights reserved.
//

#ifndef structs_h
#define structs_h

struct point {
    
    double x;
    double y;
    double z;
};

struct tensor {
    
    struct point fx;
    struct point fy;
    struct point fz;
};

#endif /* structs_h */
