function [ total_size ] = web_on(  )
%WEB_ON Summary of this function goes here
%   Detailed explanation goes here

%http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6252145

%NUMBER OF MAIN OBJECTS
main_object_number_max=212;
main_object_number=round(lognrnd(0.473844,0.68847));
%main_object_number_max_exceed_indices = main_object_number > main_object_number_max;
%main_object_number(main_object_number_max_exceed_indices) = main_object_number_max;
main_object_number=min(main_object_number,main_object_number_max);
main_object_number=max(main_object_number,1);


%MAIN OBJECT SIZE (uncompressed) (unit used: B)
main_object_size_max=8e6;                                %max: 8MB
main_object_sizes = wblrnd(28242.8,0.814944,1,main_object_number);  %unit: B
%main_object_max_exceed_indices = main_object_sizes > main_object_size_max;
%main_object_sizes(main_object_max_exceed_indices) = main_object_size_max;
main_object_sizes=min(main_object_sizes,main_object_size_max);

%NUMBER OF INLINE OBJECTS for each main object
inline_object_number_max=1920;
inline_object_number=round(exprnd(31.9291,1,main_object_number));
%inline_object_number_max_exceed_indices = inline_object_number > inline_object_number_max;
%inline_object_number(inline_object_number_max_exceed_indices) = inline_object_number_max;
inline_object_number=min(inline_object_number,inline_object_number_max);

inline_object_number_total = sum(inline_object_number);


%INLINE OBJECT SIZE (uncompressed) (unit used: B)
inline_object_size_max=8e6;                                %max: 8MB
inline_object_sizes = lognrnd(9.17979,1.24646,1,inline_object_number_total);            %unit: B
%inline_object_max_exceed_indices = inline_object_sizes > inline_object_size_max;
%inline_object_sizes(inline_object_max_exceed_indices) = inline_object_size_max;
inline_object_sizes=min(inline_object_sizes,inline_object_size_max);

total_size = sum(inline_object_sizes)+sum(main_object_sizes);
%sum(main_object_sizes')' + sum(inline_object_sizes')';

end

