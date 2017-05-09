function [isa]=logicor(a,b)

%[isa]=logicor(a,b)
%
%LOGICOR.m is a simple utility to implement a functionality in find that
%appears to be lacking. While x==1 determines where x equals 1, there
%was no easy way to do an equivalent x==[1 2 3] where you want to
%know where x equals 1, 2, or 3. While this can easily be done with a loop,
%looping is significantly slower than matrix operations in Matlab.
%Therefore, the slightly less intuitive but much faster implementation is
%enabled here. 
%
%INCLUDE:
%
%INPUTS:    A,B:        The vectors to be compared to each other. A should
%                       be a column, while b should be a row. A is also
%                       assumed to be the larger of the two. 
%
%OUTPUTS:   ISA:        The equivalent of a==b(1) or a==b(i)
%
%Created by Stephen Anthony 1/2007 U. Illinois Urbana-Champaign
%Last modified by Stephen Anthony on 09/08/2008

%Added by Stephen Anthony on 09/08/2008, can significantly improve runtime
%if there are inputs which do not match
if length(b)>1000
    b=intersect(a,b);
end

%Recoded to allow b to be a null input.
isa=false(size(a));

%Find all possibilities
for j=1:length(b)
    isa=isa | (a==b(j));
end