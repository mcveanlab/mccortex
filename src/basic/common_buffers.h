#ifndef COMMON_BUFFERS_H_
#define COMMON_BUFFERS_H_

#include "madcrowlib/madcrow_buffer.h"

madcrow_buffer(char_ptr_buf, CharPtrBuffer, char*);
madcrow_buffer(byte_buf,     ByteBuffer,    uint8_t);
madcrow_buffer(uint32_buf,   Uint32Buffer,  uint32_t);
madcrow_buffer(int32_buf,    Int32Buffer,   int32_t);
madcrow_buffer(size_buf,     SizeBuffer,    size_t);
madcrow_buffer_wipe(zsize_buf, ZeroSizeBuffer, size_t);

#endif /* COMMON_BUFFERS_H_ */
