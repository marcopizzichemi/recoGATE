#ifndef __LMFORMATS_EVENTSTORE_HPP__DEFINED__
#define __LMFORMATS_EVENTSTORE_HPP__DEFINED__

#include "formats.hpp"
#include <cstdio>

class EventStore {
public:
		EventStore();
		~EventStore();

		void addEvents(DKFZFormat *buffer, int count);
		void rewind();
		int readEvents(DKFZFormat *buffer, int count);

private:
		std::FILE *tmpFile;		
};

#endif
