function EdgeNumber = cal_edge_num(tmp,value,index,edge_seed)
n = length(value);
value = [tmp,value];
for i = 1:n
   tmp1 = (edge_seed-value(i))*(edge_seed-value(i+1));
   if tmp1<=0 && (edge_seed-value(i+1))~=0
       EdgeNumber = index(i);
       break;
   end     
end