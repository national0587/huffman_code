#include <iostream>
#include <queue>
#include <map>
#include <climits> // for CHAR_BIT
#include <iterator>
#include <algorithm>
#include <bitset>
//const int UniqueSymbols = 1 << CHAR_BIT;
const char* SampleString = "this is an example for huffman encoding";
const  int UniqueSymbols = 256;
const std::vector<int> SampleList{21,32,21,21,32,34,45,67,56,45,45,34,33};
typedef std::vector<bool> HuffCode;
typedef std::map<int, HuffCode> HuffCodeMap;

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
void bitstream::reset(){
    current_data_offset=0;
    current_bit_offset =0;
    data_size = 0;
    last_bit_offset=0;
}

void bitstream::pushbit(char c){
    if(current_bit_offset > 7){
        //次のchar配列に移動してoffsetを移動

        //std::cout << "normal flush" << std::endl;
        data.push_back(tempchar);
        ++data_size;
        tempchar = 0x0;
        current_bit_offset = 0;
        ++current_data_offset;
    }
    tempchar = (tempchar << 1) | (c == '1' ? 0x1 : 0x0);
    ++last_bit_offset;
    ++current_bit_offset;
}

void bitstream::flushbit(){
    tempchar = tempchar << (8-current_bit_offset);
    data.push_back(tempchar);
    ++data_size;
    //std::cout << "exit flush" << std::endl;
    //std::cout <<"out" << std::bitset<8>(tempchar) << std::endl;
    tempchar = 0x0;
    reset();
}

char bitstream::readbit() {
    //dataから1bitづつreturnする
    //終末が来たら\nをreturnする
    char temp;
    //data vecのsizeとdataoffsetのチェック
    if(data.size()>current_data_offset){
        //std::cout << "dataoffset "<<current_data_offset << ", bitoffset " << current_bit_offset <<std::endl;
        temp = data.at(current_data_offset);
        //bitoffset目だけを取り出す
        temp = (temp >> (7-current_bit_offset))&1;
        if(current_bit_offset == 7){
            current_bit_offset=0;
            ++current_data_offset;
        }else {
            ++current_bit_offset;
        }
    }else{
        return '\n';
    }
    return temp;
}




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

struct NodeCmp
{
    bool operator()(const INode* lhs, const INode* rhs) const { return lhs->f > rhs->f; }
};

INode* BuildTree(const int (&frequencies)[UniqueSymbols])
{
    std::priority_queue<INode*, std::vector<INode*>, NodeCmp> trees;

    for (int i = 0; i < UniqueSymbols; ++i)
    {
        if(frequencies[i] != 0)
            trees.push(new LeafNode(frequencies[i], (int)i));
    }
    while (trees.size() > 1)
    {
        INode* childR = trees.top();
        trees.pop();

        INode* childL = trees.top();
        trees.pop();

        INode* parent = new InternalNode(childR, childL);
        trees.push(parent);
    }
    return trees.top();
}





void GenerateCodes(const INode* node, const HuffCode& prefix, HuffCodeMap& outCodes)
{
    if (const LeafNode* lf = dynamic_cast<const LeafNode*>(node))
    {
        outCodes[lf->c] = prefix;
    }
    else if (const InternalNode* in = dynamic_cast<const InternalNode*>(node))
    {
        HuffCode leftPrefix = prefix;
        leftPrefix.push_back(false);
        GenerateCodes(in->left, leftPrefix, outCodes);

        HuffCode rightPrefix = prefix;
        rightPrefix.push_back(true);
        GenerateCodes(in->right, rightPrefix, outCodes);
    }
}
class HuffmanCpp
{
private:
    INode* root;
    HuffCodeMap codes;
    bitstream bs;
public:
    INode* getRoot();
    HuffCodeMap getCodeMap();
    void makeCodeBook(std::vector<int> inputVector);
    void Vec2CodeVec(std::vector<int>inputv);

};

INode *HuffmanCpp::getRoot(){
    return root;
}

HuffCodeMap HuffmanCpp::getCodeMap(){
    return codes;
}

void HuffmanCpp::makeCodeBook(std::vector<int> inputVector) {
    int frequencies[UniqueSymbols] = {0};
    for(int i=0; i<inputVector.size(); i++){
        ++frequencies[inputVector[i]];
    }
    root = BuildTree(frequencies);
    GenerateCodes(root, HuffCode(), codes);
}

void HuffmanCpp::Vec2CodeVec(std::vector<int>inputv){
    for(auto itr=inputv.begin();itr!=inputv.end();++itr){

        for(auto itr=SampleList.begin();itr!=SampleList.end();++itr){
            HuffCodeMap::iterator it = codes.find(*itr);
            if(it !=codes.end()){
                for (int i = 0; i < it->second.size(); ++i) {
//            std::cout << it->second[i] << std::endl;
                    if(it->second[i]==1){
                        bs.pushbit('1');
                    }else{
                        bs.pushbit('0');
                    }
                }
            }
        }
        bs.flushbit();
    }
}
int main()
{

    bitstream bs;

    //bs.pushbit('1');
    //終わりのおまじない


    HuffmanCpp hoge;
    HuffCodeMap codes;

    hoge.makeCodeBook(SampleList);
    codes = hoge.getCodeMap();



    for (auto itr=bs.data.begin(); itr!=bs.data.end(); ++itr) {

        std::cout << std::bitset<8>(*itr) << std::endl;

    }

    char temp;
    std::cout << "pop--->" << std::endl;
    while (1){
        temp = bs.readbit();
        if(temp=='\n') break;
        std::cout << std::hex<<(unsigned int) temp << std::endl;
    }

    //std::cout << std::bitset<8>(c)  << std::endl;

    for (HuffCodeMap::const_iterator it = codes.begin(); it != codes.end(); ++it)
    {
        std::cout << it->first << " ";
        std::copy(it->second.begin(), it->second.end(),
                  std::ostream_iterator<bool>(std::cout));
        std::cout << std::endl;
    }
    return 0;
}