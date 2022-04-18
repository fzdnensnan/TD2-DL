clc,clear all;
% network = '.\experiment\Sioux-medium-case-P=0.7.xlsx';
network = '.\experiment\Anaheim-kappa=0.1-P=0.5.xlsx';
% network = '.\CTD_EXP\test_network.xlsx';
% policyFileName = '.\CTD_EXP\policy_medium_case_18_24_0.7_10.xlsx';
policyFileName = '.\CTD_EXP\policy_Anaheim_0.1_0.5_385_48_0.7_1.xlsx';
% 1 0.001,0.001
% 2 0.01,0.001 完全收敛了
[muGraph,sigmaGraph,traversalGraph,nodeNum] = loadGraphFromFile(network);
epsilon = 0.3;
zeta = 1;
alpha = 0.001;
alphavar = 0.001;
maxEpisode = 500000;
% fullyDecisionListMeanStd = zeros(nodeNum,nodeNum);   %每一列为对应节点的decision list
% stateActionsCell = cell(1,nodeNum);
originNode = 385;
destinationNode = 48;
counter = 0;

VList = zeros(1,nodeNum);
VvarList = zeros(1,nodeNum);
[fullyDecisionListMeanStd,stateActionsCell] = initStateActionsCell(muGraph,traversalGraph,destinationNode);
plotV = [];
while counter < maxEpisode
    curNode = originNode;
    nextNode = sampleFromDLNetwork(epsilon, stateActionsCell{curNode});
    while curNode ~= destinationNode
        nextNode = sampleFromDLNetwork(epsilon, stateActionsCell{curNode});
        curReward = sampleReward(muGraph(curNode,nextNode),sigmaGraph(curNode,nextNode),'Gaussian');
        % update V
        [curUpdateV,curUpdateVvar] = Vupdate(curReward,alpha,alphavar,VList(curNode),VList(nextNode),VvarList(curNode),VvarList(nextNode));
        VList(curNode) = curUpdateV;
        VvarList(curNode) = curUpdateVvar;
        % update Q
        [curUpdateQVlist,curUpdatePolicy] = QZupdate(curReward,zeta,alpha,alphavar,stateActionsCell{curNode},VList(nextNode),VvarList(nextNode),nextNode,nodeNum);
        stateActionsCell{curNode} = curUpdateQVlist;
        fullyDecisionListMeanStd(:,curNode) = curUpdatePolicy;
        curNode = nextNode;
%         plotV = [plotV,sum(VList)];
    end
    plotV = [plotV,VList(originNode)];
    fprintf('episode: %d, V(start) = %f\n',counter, VList(originNode));
    if counter / 10000 == 0
        xlswrite(policyFileName, fullyDecisionListMeanStd);
        figure(1);
        [~,len] = size(plotV);
        x = 1 : 1 : len;
    %     plot(x, plotZ);
        set(figure(1), 'Color', 'white');% 设置图片窗口背景颜色为白色
        plot(x, plotV, 'LineWidth', 2, 'Color', [0.0118, 0.0359, 0.4824]); 
    end
    counter = counter + 1;
end
    figure(1);
    [~,len] = size(plotV);
    x = 1 : 1 : len;
%     plot(x, plotZ);
    set(figure(1), 'Color', 'white');% 设置图片窗口背景颜色为白色
    plot(x, plotV, 'LineWidth', 2, 'Color', [0.0118, 0.0359, 0.4824]); 
    

function [curUpdateQVlist,curPolicy] = QZupdate(curReward,zeta,alpha,alphavar,curQVlist,Vnext,Vvarnext,nextNode,nodeNum)
    nodeIndex = 1;
    curUpdateQVlist = curQVlist;
    [~,actionSize] = size(curQVlist);
    for i = 1 : actionSize
        if curQVlist(i).nextState == nextNode
            nodeIndex = i;
            break;
        end
    end
    delta = curReward + Vnext - curUpdateQVlist(nodeIndex).Q;
    curUpdateQVlist(nodeIndex).Q = curUpdateQVlist(nodeIndex).Q + alpha * (curReward + Vnext - curUpdateQVlist(nodeIndex).Q);
    curUpdateQVlist(nodeIndex).Qvar = curUpdateQVlist(nodeIndex).Qvar + alphavar * (delta^2 + Vvarnext - curUpdateQVlist(nodeIndex).Qvar);
    curUpdateQVlist(nodeIndex).Z = curUpdateQVlist(nodeIndex).Qvar + zeta * sqrt(curUpdateQVlist(nodeIndex).Qvar);
    
    %   rerank DL-policy
    curPolicy = zeros(nodeNum,1);
    for ei = 1 : actionSize - 1
        for ej = ei + 1 : actionSize
            if curUpdateQVlist(ei).Z > curUpdateQVlist(ej).Z
                buffer = curUpdateQVlist(ei);
                curUpdateQVlist(ei) = curUpdateQVlist(ej);
                curUpdateQVlist(ej) = buffer;
            end
        end
        curPolicy(ei,1) = curUpdateQVlist(ei).nextState;
    end
    curPolicy(actionSize,1) = curUpdateQVlist(actionSize).nextState;
    
end

function [curUpdateV,curUpdateVvar] = Vupdate(curReward,alpha,alphavar,Vcur,Vnext,Vvarcur,Vvarnext)
    curUpdateV = Vcur + alpha * (curReward + Vnext - Vcur);
    delta = curReward + Vnext - Vcur;
    curUpdateVvar = Vvarcur + alphavar * (delta^2 + Vvarnext - Vvarcur);
end

function curReward = sampleReward(mu, sigma, edgeCostDistributionModel)
    switch edgeCostDistributionModel
        case {'Gaussian','gaussian'}
            curReward = normrnd(mu,sigma);
            if curReward < 0.1
                curReward = 0.1;
            end
        otherwise
            error('no available distribution model');
     end
end

function nextNode = sampleFromDLNetwork(epsilon, curQVlist)
    nextNode = 0;
    [~,actionSize] = size(curQVlist);
    if rand() < epsilon
        randomIndex = unidrnd(actionSize);
        nextNode = curQVlist(randomIndex).nextState;
    else
        for i = 1 : actionSize
            if rand() < curQVlist(i).P
                nextNode = curQVlist(i).nextState;
                break;
            end
        end
    end
    if nextNode == 0
        nextNode = 1;
    end
end

function [muGraph,sigmaGraph,probGraph,nodeNum] = loadGraphFromFile(network)
    [xlx,~,~]=xlsread(network);
    [n,~] = size(xlx);      %
    nodeNum = xlx(n,1);
    muGraph = zeros(nodeNum,nodeNum);
    probGraph = zeros(nodeNum,nodeNum);

    % probList = zeros(n,1);
    for i = 1:n
        startState = xlx(i,1);
        endState = xlx(i,2);
        muGraph(startState,endState) = xlx(i,4);
        sigmaGraph(startState,endState) = xlx(i,5);
        probGraph(startState,endState) = xlx(i,6);
    %     probList(i) = rand();   %随机生成通行概率
    %     probGraph(Num(i,1),Num(i,2)) = probList(i);
    end
end

function [fullyDecisionListMeanStd,outputCell] = initStateActionsCell(muGraph,traversalGraph,destinationNode)
    [~,nodeNum] = size(muGraph);
    outputCell = cell(1,nodeNum);
    fullyDecisionListMeanStd = zeros(nodeNum,nodeNum);
    shortestPathList = zeros(1,nodeNum);
    c=sparse(muGraph);
    for i = 1 : nodeNum
        [dist,~,~] = graphshortestpath(c,i,destinationNode);
        shortestPathList(i) = dist;
    end
    
    for i = 1 : nodeNum
%         VfuncList(i) = 0;
        if i ~=destinationNode && shortestPathList(i) ~= inf
            bufQVlist = [];        %Q(s,a)、J(s,a)、
            for j = 1 : nodeNum
               if  muGraph(i,j)~=0 && shortestPathList(j) ~= inf
%                    bufQPJ.Q = 0;
                   bufQV.Q = shortestPathList(j) + muGraph(i,j);
                   bufQV.Qvar = 0;
%                    bufQV.V = 0;
%                    bufQV.Vvar = 0;
                   bufQV.P = traversalGraph(i,j);
                   bufQV.nextState = j;
                   bufQV.Z = bufQV.Q;
                   bufQV.preZ = 0;
                   bufQVlist = [bufQVlist,bufQV];
               end
            end
            [~,sizeBufList] = size(bufQVlist);
            if sizeBufList > 0
                for ei = 1 : sizeBufList - 1
                    for ej = ei + 1 : sizeBufList
                        if bufQVlist(ei).Q > bufQVlist(ej).Q
                            buffer = bufQVlist(ei);
                            bufQVlist(ei) = bufQVlist(ej);
                            bufQVlist(ej) = buffer;
                        end
                    end
                    fullyDecisionListMeanStd(ei,i) = bufQVlist(ei).nextState;
                end
                 fullyDecisionListMeanStd(sizeBufList,i) = bufQVlist(sizeBufList).nextState;
            end
            outputCell{i} = bufQVlist;
        else
            bufQV.Q = 0;
            bufQV.Qvar = 0;
%             bufQV.V = 0;
%             bufQV.Vvar = 0;
            bufQV.P = 0;
%             bufQV.nextState = j;
            bufQV.Z = bufQV.Q;
            bufQV.preZ = 0;
            outputCell{i} = bufQV;
        end
    end
end
