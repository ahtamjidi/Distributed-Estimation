function calc_cent_update_gold()
run_cent_filter()
end
function run_cent_filter()
global opt_dist

% initialization
x_bar = opt_dist.result.initial.x_bar(:,1);
P_bar =  opt_dist.result.initial.P_bar(:,:,1);

x_update = x_bar;
P_update = P_bar;

for i_step = 1 : opt_dist.i_step
    % predict
    [x_pred,P_pred] = f(x_update,P_update,0);
    Y_bar = pinv(P_pred);
    y_bar = Y_bar*x_pred;
    % calc update
    delta_I = zeros(size(Y_bar));
    delta_i = zeros(size(y_bar));
    
    %     [size_group,nComponents,members] = networkComponents_gold(opt_dist.Graph_History{i_step});
    %     neigbours=members{nComponents(i_agent)};
    [delta_I,delta_i] = collect_observations_cent( i_step,delta_I,delta_i);
    
    %     for i_neighbour = 1:numel(neigbours)
    %         H = opt_dist.result.obs.H{i_step,neigbours(i_neighbour)};
    %         z = opt_dist.sim.obs.z{i_step,neigbours(i_neighbour)};
    %         delta_I = delta_I + 1/(opt_dist.sim.obs.r_var{i_step,neigbours(i_neighbour)})*(H'*H);
    %         delta_i = delta_i +  1/(opt_dist.sim.obs.r_var{i_step,neigbours(i_neighbour)})*H'*z;
    %     end
    % update
    Y_update = Y_bar + delta_I;
    y_update = y_bar + delta_i;
    
    P_update = pinv(Y_update);
    x_update = P_update*y_update;
    
    opt_dist.result.estcent{i_step}.Y_bar = Y_update;
    opt_dist.result.estcent{i_step}.y_bar = y_update;
    opt_dist.result.estcent{i_step}.P_bar = P_update;
    opt_dist.result.estcent{i_step}.x_bar = x_update;
    
    
end


x_update
P_update

% opt_dist.connection_history = connection_history;

end
function [delta_I,delta_i] = collect_observations_cent(i_step,delta_I,delta_i)

global  opt_dist
for i_agent = 1: opt_dist.nAgents
    H = opt_dist.result.obs.H{i_step,i_agent};
    z = opt_dist.sim.obs.z{i_step,i_agent};
    delta_I = delta_I + 1/(opt_dist.sim.obs.r_var{i_step,i_agent})*(H'*H);
    delta_i = delta_i +  1/(opt_dist.sim.obs.r_var{i_step,i_agent})*H'*z;
end
end