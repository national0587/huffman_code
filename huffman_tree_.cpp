
#include "huffman_tree_.h"
//const int UniqueSymbols = 1 << CHAR_BIT;

using namespace huffman_test;
const char* SampleString = "this is an example for huffman encoding";
const  int UniqueSymbols = 256;

void bitstream::reset(){
    current_data_offset=0;
    current_bit_offset =0;
    //data_size = 0;
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
//    std::cout << "--readbit finish--" << std::endl;
    return temp;
}




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



INode *HuffmanCpp::getRoot(){
    return root;
}

HuffCodeMap HuffmanCpp::getCodeMap(){
    return codes;
}


void HuffmanCpp::makeCodeBook(std::vector<int> &inputVector) {
    int frequencies[UniqueSymbols] = {0};
    for(int i=0; i<inputVector.size(); i++){
        ++frequencies[inputVector[i]];
    }
    root = BuildTree(frequencies);
    GenerateCodes(root, HuffCode(), codes);
}

void HuffmanCpp::Vec2CodeVec(std::vector<int> &inputv){
    //vectorをbitstreamに変換するメソッド
    //入力と出力は参照渡し

    for(auto itr=inputv.begin();itr!=inputv.end();++itr){
        HuffCodeMap::iterator it = codes.find(*itr);
        if(it !=codes.end()){
            for (int i = 0; i < it->second.size(); ++i) {
            std::cout << it->second[i] << std::endl;
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
void HuffmanCpp::FollowTree(std::vector<int> &v) {
    //符号をもとにtreeを辿る
    char chartemp;
    INode *nodetemp;
    nodetemp=root;
    while(1){
        chartemp = bs.readbit();

//        std::cout << std::hex <<(unsigned int)chartemp << std::endl;
        if(chartemp=='\n'){std::cout << "break"<< std::endl; break;}
        else{
//            std::cout << "oooou" << std::endl;
            if (const LeafNode* lf = dynamic_cast<const LeafNode*>(nodetemp))
            {
//                std::cout << "o1u" << std::endl;
            }
            else if (const InternalNode* in = dynamic_cast<const InternalNode*>(nodetemp))
            {
//                std::cout << "ou" << std::endl;
                if(chartemp==0x0){

                }else{

                }
            }
        }

    }
}

bitstream HuffmanCpp::getBitstream() {
    return bs;
}


int main()
{
    std::vector<int> SampleList{21,32,21,21,32,34,45,67,56,45,45,34,33};
    bitstream bs;

    //bs.pushbit('1');
    //終わりのおまじない

    HuffmanCpp hoge;
    HuffCodeMap codes;

    hoge.makeCodeBook(SampleList);
    codes = hoge.getCodeMap();

    hoge.Vec2CodeVec(SampleList);

    std::cout << "hoge" << std::endl;
    std::vector<int> output;
    hoge.FollowTree(output);
    std::cout << "end" << std::endl;
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