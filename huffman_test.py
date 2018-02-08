#!/usr/bin/env python
# -*- coding: utf-8 -*-
import huffman_tree as hu

def main():
    f = hu.PyHuffman()
    l = [1,2,3,4,5,5,6,4,3,2,1,1,2,3,4,5,3,2,3,4,5,2,2,2]
    f.makeCodeBook(l)
    print f.V2CV(l)


if __name__ == '__main__':
    main()