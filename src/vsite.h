//
// Created by Zhang, Xiaohua on 5/5/22.
//

#ifndef DDCMD_VSITE_H
#define DDCMD_VSITE_H

void vsite_Update(GROUP *g, int mode, STATE *state, double time_in, double dt_in);
void vsite_velocityUpdate(int mode, int k,GROUP *g, STATE *state, double time, double dt);

void vsite_construction(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e);
void vsite_force(SYSTEM*sys, CHARMMPOT_PARMS *parms, ETYPE *e);

/*
void coor_vsite1(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2);
void coor_vsite2(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3);
void coor_vsite2fd(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3);
void coor_vsite3(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4);
void coor_vsite3out(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4);

void force_vsite1(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2);
void force_vsite3(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4);
void force_vsite3out(STATE* state, VSITE_CONN *vsiteConn, int index1, int index2, int index3, int index4);
*/

#endif //DDCMD_VSITE_H
