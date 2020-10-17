#pragma once

struct edgepoint
{
	int Z;
	int LR; //1 means z + 1 is empty 2 means z -1 is empty 0 means both
};

class LLi
{
public:
	void add(int i);
	bool search(int i);
	~LLi();
	

private:
	struct node
	{
		int data;
		node* next;
	};
	node* head = 0;
};

class LLE
{
public:
	void add(edgepoint i);
	bool search(int i);
	edgepoint* dump(int *size);
	~LLE();

private:
	struct node
	{
		edgepoint data;
		node* next;
	};
	node* head = 0;
};