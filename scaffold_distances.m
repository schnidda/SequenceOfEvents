%% 180613 - using dijkstra algorithm offered from matlab to calculate the mean base distance for the folding

%function abc = abctest(x)
%abc=x
%end

function [radiusofgyration]=scaffold_distances(scaffold_length, termini);%, A, C);%gyradius(scaffold_length, termini)
% scaffold_length tells the algorithm the size of the arrays, termini
% should include the bases which will be connected and the time at which
% this loop is closed
%e.g.
% from_base to_base mean_time
%       456     3820    35.2
%% first create two arrays, scaffold: connections on the scaffold, distance: distance between bases
%then sort the termini according to time

loops=sortrows(termini,3);

scaffold=zeros(scaffold_length,scaffold_length);

distance=zeros(scaffold_length,scaffold_length);


assignin('base', 'scaffold', scaffold);
assignin('base', 'distance', distance);

%apply first entries

for i = 1:scaffold_length    
    if i~=7560
    scaffold(i,i+1)=1;
    distance(i,i+1)=1;
    end
    
    if i~=1
    scaffold(i,i-1)=1;
    distance(i,i-1)=1;
    end       
end


%%looping through the incorporation times in order to change the distance
%%map and calculate radius of gyration

radiusofgyrationintime=zeros(size(loops,1));

for i = 1:size(loops,1)
   scaffold(loops(i,1), loops(i,2)) = 1;%first argument should be i
   distance(loops(i,1), loops(i,2)) = 1;%first argument should be i
   
   %[costs,paths] = dijkstra(scaffold,distance);
   %%%%%%%
   [costs] = dijkstramod(scaffold,distance);
   
   %costs2=dijkstramod(A,C);
   %assignin('base', 'costs2', costs2);
   

   costs4=costs.^2;%costs2.^2;

   radiusofgyration = sqrt(1/2/(scaffold_length)^2*sum(sum(costs4)))
   
   radiusofgyrationintime(i)=radiusofgyration;
end


assignin('base', 'radiusofgyrationintime', radiusofgyrationintime);

end
