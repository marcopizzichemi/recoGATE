#ifndef __ELM2_HPP__
#define __ELM2_HPP__

#include <sys/types.h>

namespace ELM2 {
	struct EventFormat {
		double ts;
		u_int8_t random;
		float d;
		float yozRot;		
		float x1;
		float y1;
		float z1;
		float e1;
	        u_int8_t n1;
		float x2;
		float y2;
		float z2;
		float e2;
	        u_int8_t n2;
		float dt;
	} __attribute__((__packed__));
}
#endif
