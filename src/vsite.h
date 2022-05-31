//
// Created by Zhang, Xiaohua on 5/5/22.
//

#ifndef DDCMD_VSITE_H
#define DDCMD_VSITE_H

void vsite_Update(GROUP *g, int mode, STATE *state, double time_in, double dt_in);
void vsite_velocityUpdate(int mode, int k,GROUP *g, STATE *state, double time, double dt);

#endif //DDCMD_VSITE_H
