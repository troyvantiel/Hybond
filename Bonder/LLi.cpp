#include "stdafx.h"
#include "LLI.h"


void LLi::add(int i)
{
	node* toAdd = new node;
	(*toAdd).data = i;
	(*toAdd).next = head;
	head = toAdd;
}

bool LLi::search(int i)
{
	node *next = head;
	while (next != 0)
	{
		if ((*next).data == i)
			return true;
		next = (*next).next;
	}
	return false;
}

LLi::~LLi()
{
	node *current = head, *lastUsed;
	while (current)
	{
		lastUsed = current;
		current = (*current).next;
		delete lastUsed;
	}
}

void LLE::add(edgepoint i)
{
	node* toAdd = new node;
	(*toAdd).data = i;
	(*toAdd).next = head;
	head = toAdd;
}

bool LLE::search(int i)
{
	node *next = head;
	while (next != 0)
	{
		if ((*next).data.Z == i)
			return true;
		next = (*next).next;
	}
	return false;
}

edgepoint* LLE::dump(int* size)
{
	*size = 0;
	node *next = head;
	while (next != 0)
	{
		(*size)++;
		next = (*next).next;
	}
	if (*size == 0)
		return 0;
	edgepoint* result = new edgepoint[(*size)];
	next = head;
	for (size_t i = 0; i < (*size); i++)
	{
		result[i] = (*next).data;
		next = (*next).next;
	}
	return result;
}

LLE::~LLE()
{
	node *current = head, *lastUsed;
	while (current)
	{
		lastUsed = current;
		current = (*current).next;
		delete lastUsed;
	}
}
