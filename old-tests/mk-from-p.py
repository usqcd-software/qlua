#!/usr/bin/env python
import tables
import numpy

foo = tables.openFile("from-p.h5", "w")
x = numpy.array([1+2j, 4+3j, 6.25-9.5j])
foo.createArray("/", "c1", x)
foo.create_group("/", "subdirA")
foo.createArray("/subdirA", "r1", numpy.array([1,2,6,5]))
foo.createArray("/subdirA", "r2", numpy.array([[-1.126,5.375], [1,2], [3.5,7]]))
foo.createArray("/subdirA", "c2", numpy.array([1+2j,6j,-5]))
foo.createArray("/", "s1", "The quick brown fox") ## does not match qlua notion of strings

foo.close()
