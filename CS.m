classdef CS < ALGORITHM
    % <multi> <large> <real>
    % conversed sequences based framework
    % alpha             ---  0.5  ---
    % rate              ---  0.8  ---
    % window            ---  175  ---

    methods
        function main(Algorithm, Problem)
            % MOEAD/DE的参数
            delta = 0.9; nr = 2;

            % 算法参数设置
            [alpha,rate,cur_window_change] = Algorithm.ParameterSet(0.5,0.8,175);

            % 算法通用设置 如初始化 距离向量 邻居向量等
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);
            Population = Problem.Initialization();
            Z = min(Population.objs,[],1);

            directions = zeros(Problem.N,Problem.D);

            cnt = 0;
            T = [];
            while Algorithm.NotTerminated(Population)

                cur_alpha = alpha*exp(-((Algorithm.pro.FE.^2*log(alpha*10000)/(rate^2+Algorithm.pro.maxFE^2))));

                bak_Population = Population;
                for c = 1 : Problem.N
                    if rand < delta
                        P = B(c,randperm(end));
                    else
                        P = randperm(Problem.N);
                    end
                    offspring = OperatorDE(Population(c),Population(P(1)),Population(P(2)));
                    Z = min(Z,offspring.obj);
                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                    g_new = max(repmat(abs(offspring.obj-Z),length(P),1).*W(P,:),[],2);
                    Population(P(find(g_old>=g_new,nr))) = offspring;
                end



                if ~all(directions==0) & Algorithm.pro.FE < rate*Algorithm.pro.maxFE
                    cnt = cnt + 1;

                    tic;
                    new_directions = Population.decs()-bak_Population.decs();
                    % 找到两次更新的方向不一致的
                    idx1 = directions>0;
                    idx2 = new_directions<0;
                    idx = idx1&idx2;

                    idx1 = directions<0;
                    idx2 = new_directions>0;
                    idx = idx|(idx1&idx2);
                    [~,composition_idx] = sort(abs(Population.decs()-bak_Population.decs()),2);
                    % 对齐窗口大小
                    for i =1:Problem.N
                        if sum(idx(i,:)) < cur_window_change
                            j = 1;
                            while true
                                if idx(i,composition_idx(i,j)) == 0
                                    idx(i,composition_idx(i,j)) = 1;
                                end
                                j=j+1;
                                if sum(idx(i,:)) == cur_window_change
                                    break;
                                end
                            end
                        end
                        if sum(idx(i,:)) > cur_window_change
                            j = Problem.D;
                            while true
                                if idx(i,composition_idx(i,j)) == 1
                                    idx(i,composition_idx(i,j)) = 0;
                                end
                                j=j-1;
                                if sum(idx(i,:)) == cur_window_change
                                    break;
                                end
                            end
                        end
                    end
                    decs = Population.decs;
                    for i=1:Problem.N

                        decs(i,~idx(i,:)) = decs(i,~idx(i,:))+cur_alpha.*directions(i,~idx(i,:));

                        if sum(idx(i,:)) ~= 0
                            if rand < delta
                                P = B(i,randperm(end));
                            else
                                P = randperm(Problem.N);
                            end
                            copy_d = Problem.D; copy_lower = Problem.lower;copy_upper = Problem.upper;
                            Problem.D = cur_window_change;Problem.upper = Problem.upper(idx(i,:));Problem.lower = Problem.lower(idx(i,:));
                            cdoffspring = OperatorDE(decs(i,idx(i,:)),decs(P(1),idx(P(1),:)),decs(P(2),idx(P(2),:)));



                            %                                 sim_parental_idx = [i, P(1), P(2)];
                            %                                 sim_parental = [decs(i,idx(i,:));decs(P(1),idx(P(1),:));decs(P(2),idx(P(2),:))];
                            %                                 distance = sum((cdoffspring - sim_parental).^2,2);
                            %                                 [~, parent_idx] = min(distance);
                            %                                 parent = sim_parental_idx(parent_idx);
                            decs(i,idx(i,:)) = cdoffspring;
                                                            Problem.D=copy_d;Problem.upper=copy_upper;Problem.lower=copy_lower;
                        end


                        offspring = SOLUTION(decs(i,:));
                        Z = min(Z,offspring.objs);
                        g_old = max(abs(Population(i).objs - Z).*W(i,:),[],2);
                        g_new = max(abs(offspring.objs - Z).*W(i,:),[],2);
                        if g_old  > g_new
                            Population(i) = offspring;
                        end
                    end
                    T(cnt) = toc;
                end
                directions = directions*0+Population.decs() - bak_Population.decs();


                %             if Problem.FE >= Problem.maxFE
                %                 al = 0.5;
                %                 be = 0.25;
                %                 disp("alg finished")
                %                 Tt = sum(T);
                %                 tao = Tt/Algorithm.metric.runtime;
                %                 Fe = (Problem.FE - Problem.maxFE)/Problem.maxFE;
                %                 P = IGD(Population,Problem.optimum);
                %                 f = al*P + be*tao + (1-al-be)*Fe;
                %
                %                 res = [Problem.name,'\t',num2str(f),'\t',num2str(tao),'\n'];
                %                 fid = fopen('CS.txt','a+');
                %                 fprintf(fid,'%s',res);
                %                 fclose(fid);
                %             end

            end


        end
    end
end

% 0.5*exp(-((x.^2*log(5000)/(0.8^2+100000^2))))
