#include "crc32.h"

#include <stdio.h>

static unsigned char crc32_lut_inited = 0;
static unsigned  crc32_lut[256];

static unsigned reflect(unsigned long crc, int bitnum)
{
   // reflects the lower 'bitnum' bits of 'crc'
   unsigned i, j=1, crcout=0;

   for (i=(unsigned)1<<(bitnum-1); i; i>>=1)
   {
      if (crc & i) crcout|=j;
      j<<= 1;
   }
   return (crcout);
}

static void init_crc32_table(void)
{
   if(crc32_lut_inited) return;

   int i, j;
   unsigned bit, crc;
   for (i=0; i<256; i++)
   {
      crc=(unsigned )i;
      crc=reflect(crc, 8);
      crc<<= 32-8;

      for (j=0; j<8; j++)
      {
         bit = crc & 0x80000000;
         crc<<= 1;
         if (bit) crc^= 0x4c11db7;
      }

      crc = reflect(crc, 32);
      crc32_lut[i]= crc;
   }
   crc32_lut_inited = 1;
}

unsigned checksum_crc32(unsigned char * data, unsigned len)
{
   unsigned i, j, c, bit;
   unsigned crc = 0xffffffff;
   for (i=0; i<len; i++)
   {
      c = (unsigned)*data++;
      c = reflect(c, 8);

      for (j=0x80; j; j>>=1)
      {
         bit = crc & 0x80000000;
         crc<<= 1;
         if (c & j) bit^= 0x80000000;
         if (bit) crc^= 0x4c11db7;
      }
   }
   crc=reflect(crc, 32);
   crc^= 0xffffffff;
   return crc ;
}



unsigned checksum_crc32_table(unsigned char * data, unsigned len)
{
   if(!crc32_lut_inited) {init_crc32_table();}

   unsigned crc = 0xffffffff;
   /*crc = reflect(crc, 32);*/
   while(len--)
   {
      crc = (crc >> 8) ^ crc32_lut[ (crc & 0xff) ^ *data++];
   }
   crc^= 0xffffffff;
   return(crc);
}


/*

int main(int argc, char ** argv)
{
   char * test = "123456789";
   unsigned  crc32;

   crc32  = checksum_crc32(test,strlen(test));
   printf("Bit By Bit:\n");
   printf("\tHex: %u\n",crc32);
   printf("\tInt: %x\n",crc32);

   crc32  = checksum_crc32_table(test,strlen(test));
   printf("Table:\n");
   printf("\tHex: %u\n",crc32);
   printf("\tInt: %x\n",crc32);
}





*/


/* Local Variables: */
/* tab-width: 3 */
/* End: */
