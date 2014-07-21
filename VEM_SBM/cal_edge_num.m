function EdgeNumber = cal_edge_num(tmp,pro,x,edge_seed)
n = length(pro);
pro = [tmp,pro];
x = [0,x];
EdgeNumber = 0;
for i = 1:n
   tmp1 = (edge_seed-pro(i))*(edge_seed-pro(i+1));
   if tmp1<=0 && (edge_seed-pro(i+1))~=0
       EdgeNumber = x(i);
       break;
   end
end