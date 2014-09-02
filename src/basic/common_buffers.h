#ifndef COMMON_BUFFERS_H_
#define COMMON_BUFFERS_H_

#include "objbuf_macro.h"

create_objbuf(char_ptr_buf, CharPtrBuffer, char*);
create_objbuf(byte_buf,     ByteBuffer,    uint8_t);
create_objbuf(uint32_buf,   Uint32Buffer,  uint32_t);
create_objbuf(int32_buf,    Int32Buffer,   int32_t);
create_objbuf(size_buf,     SizeBuffer,    size_t);

#endif /* COMMON_BUFFERS_H_ */
