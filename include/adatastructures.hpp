#ifndef _GLIBCXX_IOSTREAM
#include <iostream>
#endif

#define __ADATASTRUCTURES_H

using namespace std;

template <class T>
struct Node
{
	T data;
	Node* next;
};

template <class T>
class CircList
{
private:
	Node<T>* start;
	int size;

public:
	CircList(T n)
	{
		start = new Node<T>;
		start->data = n;
		start->next = start;
		size = 1;
	}

	//Get the node which is "n" nodes away from start
	Node<T>* traverse(int n)
	{
		if(n < 0) { cout << "! CircList: traverse: invalid argument!\n"; return nullptr;}
		Node<T>* cur;
		cur = start;
		int i = 0;

		while(true)
		{
			if(i == n) break;
			cur = cur->next;
			i++;
		}

		return cur;
	}

	// Get last node (before start)
	Node<T>* last()
	{
		Node<T>* cur;
		cur = start;

		while(cur->next != start)
			cur = cur->next;

		return cur;
	}

	void push(T x)
	{
		size++;
		Node<T>* nnew = new Node<T>;
		nnew->data = x;
		nnew->next = start;
		last()->next = nnew;
	}

	// return pointer to first node containing data x
	Node<T>* find(T x)
	{
		if (start->data == x) return start;
		Node<T>* cur = start->next;
		while(cur != start)
		{
			if(cur->data == x) return cur;
			cur = cur->next;
		}
		cout << "! CircList: find: Data not found in list!\n";
		return nullptr;
	}

	void print_list()
	{
		Node<T>* cur;
		cur = start;
		cout << start->data << " ";

		while(cur->next != start)
		{	cur = cur->next;
			cout << cur->data << " ";
		}
		cout << endl;
	}

	~CircList()
	{
		Node<T> *cur, *nnext;
		cur = start->next;
		while(cur != start)
		{
			nnext = cur->next;
			delete cur;
			cur = nnext;
		}
		delete start;
	}
};

// This function is a cyclic permutation of consecutive integers from 'start' to 'end' (inclusive). It returns the integer (between 'start' and 'end') that is 'off' integers away from 'n' in the cyclic order.
int perm(int start, int end, int n, int off)
{
	if(n > end) { cout << "Permutation point error!\n"; return 0; }
	if(off == 0) return n;

	CircList<int> list(start);
	for(int i = start+1; i <= end; i++)
		list.push(i);

	Node<int>* nn = list.find(n);
	Node<int>* cur = nn;
	for(int i = 0; i < off; i++)
		cur = cur->next;
	return cur->data;
}

template <class T>
class MergeSort
{
	T** arrs;
	int length;
	int num_arr;
	int nsort;
public:
	MergeSort(T** arrays, int length_arrays, int num_arrays, int index_of_sorting_array)
		: arrs(arrays), length(length_arrays), num_arr(num_arrays), nsort(index_of_sorting_array)
	{  }

	void merge(T** a, T** b, T** fin, int len)
	{
		int i = 0, j = 0;
	}

	void mergesort()
	{
	}
};
