//
// Created by hiroo on 18/02/02.
//

#ifndef HUFFMAN_CYTHON_HUFFMAN_TREE_H
#define HUFFMAN_CYTHON_HUFFMAN_TREE_H

#include <iostream>
#include <queue>
#include <map>
//#include <climits> // for CHAR_BIT
#include <iterator>
#include <algorithm>
#include <bitset>

namespace huffman_test{

    typedef std::vector<bool> HuffCode;
    typedef std::map<int, HuffCode> HuffCodeMap;

class INode
{
public:
    const int f;

    virtual ~INode() {}

protected:
    INode(int f) : f(f) {}
};

class InternalNode : public INode
{
public:
    INode *const left;
    INode *const right;

    InternalNode(INode* c0, INode* c1) : INode(c0->f + c1->f), left(c0), right(c1) {}
    ~InternalNode()
    {
        delete left;
        delete right;
    }
};

class LeafNode : public INode
{
public:
    const int c;

    LeafNode(int f, int c) : INode(f), c(c) {}
};

class bitstream{
private:
    char tempchar=0x0;
    int current_data_offset;
    int current_bit_offset;
public:
    //reset
    bitstream():
            current_data_offset(0)
            ,current_bit_offset(0)
            ,data_size(0)
            ,last_bit_offset(0)
    {}
    void reset();
    void pushbit(char c);
    void flushbit();
    char readbit();
    std::vector<char> data;
    int data_size; //sizeof data
    int last_bit_offset; //size of bit
};

class HuffmanCpp
{
private:
    INode* root;
    HuffCodeMap codes;
    bitstream bs;
public:
    HuffmanCpp(){};
    ~HuffmanCpp(){
    };
    INode* getRoot();
    HuffCodeMap getCodeMap();
    bitstream getBitstream();
    void makeCodeBook(std::vector<int> &inputVector);
    void Vec2CodeVec(std::vector<int> &inputv);
    void FollowTree(std::vector<int> &);

};
}
#endif //HUFFMAN_CYTHON_HUFFMAN_TREE_H
