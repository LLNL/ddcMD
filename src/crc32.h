#ifndef CRC32_H
#define CRC32_H

#ifdef __cplusplus
extern "C"{
#endif

unsigned checksum_crc32(unsigned char * data, unsigned len);
unsigned checksum_crc32_table(unsigned char * data, unsigned len);

#ifdef __cplusplus
}
#endif

#endif // #ifndef CRC32_H


/* Local Variables: */
/* tab-width: 3 */
/* End: */
