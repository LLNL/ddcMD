/* 
 * File:   bioGid.h
 * Author: zhang
 *
 * Created on May 28, 2015, 2:23 PM
 */

#ifndef BIOGID_H
#define	BIOGID_H

#include "gid.h"

    // gid type MOL              RESID   GROUP  ATOM
    //          |------32------|---16---|--8--|--8--|
static const int  molShift=32; 
static const gid_type atmMask=255ull;
static const gid_type atmgrpMask=65535ull;
static const gid_type grpMask=(255ull<<8);
static const gid_type resMask=(65535ull<<16);
static const gid_type molMask=(4294967295ull<<32);
static const gid_type molResMask = 0xffffffffffff0000;
static const gid_type totMask=(4294967295ull<<32) | (65535ull<<16) | 65535ull;
static const gid_type molresMask=~(65535ull);
#endif	/* BIOGID_H */

