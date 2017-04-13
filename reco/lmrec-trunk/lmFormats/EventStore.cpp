#include "EventStore.hpp"
#include <cstdio>
#include <assert.h>

EventStore::EventStore()
{
	tmpFile = std::tmpfile();
	assert(tmpFile != NULL);
}

EventStore::~EventStore()
{
	fclose(tmpFile);
}

void EventStore::addEvents(DKFZFormat *buffer, int count)
{
	int r = std::fwrite((void *)buffer, sizeof(DKFZFormat), count, tmpFile);
	assert (r == count);
}

void EventStore::rewind()
{
	std::rewind(tmpFile);
}

int EventStore::readEvents(DKFZFormat *buffer, int count)
{
	int r = std::fread((void *)buffer, sizeof(DKFZFormat), count, tmpFile);
	return r;
}
