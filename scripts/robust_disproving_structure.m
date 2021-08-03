%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symbolic Parameterized Unitary Residual Calculator
% August 3rd, 2021
% Max Aksel Bowman (mbowman@anl.gov)
% 
% This code is designed to return a vector of complex-valued symbolic
% residuals representing the element-wise differences between elements of a
% unitary matrix representing a quantum circuit constructed in a
% layer-based fashion (shown below) parameterized by 3n+3 real parameters
% (where n is the number of multi-qubit operations required to represent
% the circuit) and some target unitary matrix. A layer-based circuit may
% look as follows:
%       ------           ------                   ------
% q0 ---| U3 | --- * --- | U3 | --- * --- ... --- | U3 |
%       ------     |     ------     |             ------
%       ------     |     ------     |             ------
% q1 ---| U3 | --- x --- | U3 | --- | --- ... --- | U3 |
%       ------     |     ------     |             ------
%       ------     |     ------     |             ------
% q2 ---| U3 | --- x --- | U3 | --- x --- ... --- | U3 |
%       ------           ------                   ------
% Qubits are indexed 0 through 2. Each layer consists of three U3 gates
% on each qubit followed by one of three multi-qubit gates: a CNOT between
% qubits 0 and 1, a CNOT between qubits 0 and 2, and a quantum FAN-OUT
% operation controlled by qubit 0. This restricted set of multi-qubit
% operations is required to respect a linear qubit topology.
%
% Operations are defined using Little Endian.
%
% The user may alter the target unitary and the layer_instructions array.
% To construct a 3-layer circuit with the all FAN-OUT operations, for
% example, the layer_instructions array should be set to [2 2 2].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global residuals;
main();

function main
    % Program entry
    global residuals;
    tic;

    % Define symbolic variables and multi-qubit operations
    layer_instructions = [2, 1, 0];
    thetas  = sym('theta',  [1, 3*(numel(layer_instructions)+1)]);
    phis    = sym('phi',    [1, 3*(numel(layer_instructions)+1)]);
    lambdas = sym('lambda', [1, 3*(numel(layer_instructions)+1)]);
    
    % Define target unitary matrix (may be complex)
    target_unitary = [1 0 0 0 0 0 0 0;
                      0 1 0 0 0 0 0 0;
                      0 0 1 0 0 0 0 0;
                      0 0 0 1 0 0 0 0;
                      0 0 0 0 1 0 0 0;
                      0 0 0 0 0 1 0 0;
                      0 0 0 0 0 0 0 1;
                      0 0 0 0 0 0 1 0];
                  
    % Construct unitary from symbolic variables and multi-qubit instructions
    actual_unitary = eye(8);
    
    % Search for previously computed unitarys that could speed up
    % computation
    prev_sol_ind = numel(layer_instructions);
    while prev_sol_ind >= 1
        fname = strcat("U", num2str(layer_instructions(1:prev_sol_ind), "%d"), ".mat");
        if exist(fname, "file") == 2
            disp("Found file...");
            disp(fname);
            disp("Loading precomputed unitary...");
            load(fname, "actual_unitary");  % will contain actual_unitary
            break;
        end
        prev_sol_ind = prev_sol_ind-1;
    end
    
    if prev_sol_ind == 0
        % No prior basis found
        % Add a layer of single-qubit U3 gates
        for qubit=0:2
            theta  = thetas(qubit+1);
            phi    = phis(qubit+1);
            lambda = lambdas(qubit+1);
            actual_unitary = actual_unitary * u3_gate(qubit, theta, phi, lambda);
        end
    end
    
    for layer_index=prev_sol_ind+1:numel(layer_instructions)
        % Add a multi-qubit gate
        actual_unitary = actual_unitary * linear_mq_gate(layer_instructions(layer_index));

        % Add a layer of single-qubit U3 gates
        for qubit=0:2
            theta  = thetas(layer_index*3+1+qubit);
            phi    = phis(layer_index*3+1+qubit);
            lambda = lambdas(layer_index*3+1+qubit);
            actual_unitary = actual_unitary * u3_gate(qubit, theta, phi, lambda);
        end
        
        % Save unitary at current layer for future reference (DP)
        save(strcat("U", num2str(layer_instructions(1:layer_index), "%d"), ".mat"), "actual_unitary");
    end
    
    % Calculate and vectorize element-wise differences between target and
    % symbolic unitary matrices
    difference_matrix = target_unitary - actual_unitary;
    residuals = difference_matrix(:);
    save(strcat('residuals_', num2str(layer_instructions, "%d"),'.mat'), 'residuals');

    % Print residuals and (maybe) their gradients to the console
    residuals_file = fopen(strcat("residuals_", num2str(layer_instructions, "%d"), ".txt"), 'w');
    for ri=1:numel(residuals)
        fprintf(residuals_file, "%s\n", residuals(ri));
        % sprintf("%dth residual: %s\n", ri, residuals(ri));
        % sprintf("%dth residual gradient: ", ri)
        % disp(gradient(residuals(ri)));
    end
    fclose(residuals_file);
    
    toc;
end

function u3_mat = u3_gate(qubit, theta, phi, lambda)
% u3_gate Returns the 8x8 unitary matrix of a single-qubit gate
% Inputs:
%    qubit: an integer of value 0, 1, or 2
%    theta: a floating-point value in the range [0, 2pi)
%    phi: a floating-point value in the range [0, 2pi)
%    lambda: a floating-point value in the range [0, 2pi)
%
% Outputs:
%    mq_mat: an 8x8 unitary matrix

    u3_mat = [cos(theta/2)             -exp(1j*lambda)*sin(theta/2);
              exp(1j*phi)*sin(theta/2)  exp(1j*(phi+lambda))*cos(theta/2)];

    if qubit == 0
        u3_mat = kron(u3_mat, kron(eye(2), eye(2)));
    elseif qubit == 1
        u3_mat = kron(eye(2), kron(u3_mat, eye(2)));
    elseif qubit == 2
        u3_mat = kron(eye(2), kron(eye(2), u3_mat));
    else
        return;
    end
end

function mq_mat = linear_mq_gate(gate_type)
% linear_mq_gate Returns the 8x8 unitary matrix of a multi-qubit gate
% Inputs:
%    gate_type: an integer of value 0, 1, or 2
%
% Outputs:
%    mq_mat: an 8x8 unitary matrix

    if gate_type == 0
        % Apply a CNOT gate between qubits 0 and 1 on a three-qubit circuit
        mq_mat = [1 0 0 0 0 0 0 0;
                  0 1 0 0 0 0 0 0;
                  0 0 1 0 0 0 0 0;
                  0 0 0 1 0 0 0 0;
                  0 0 0 0 0 0 1 0;
                  0 0 0 0 0 0 0 1;
                  0 0 0 0 1 0 0 0;
                  0 0 0 0 0 1 0 0];
        return;
    elseif gate_type == 1
        % Apply a CNOT gate between qubits 0 and 2 on a three-qubit circuit
        mq_mat = [1 0 0 0 0 0 0 0;
                  0 1 0 0 0 0 0 0;
                  0 0 1 0 0 0 0 0;
                  0 0 0 1 0 0 0 0;
                  0 0 0 0 0 1 0 0;
                  0 0 0 0 1 0 0 0;
                  0 0 0 0 0 0 0 1;
                  0 0 0 0 0 0 1 0];
         return;
    elseif gate_type == 2
        % Apply a FAN-OUT gate with qubit 0 as the control qubit
        mq_mat = [1 0 0 0 0 0 0 0;
                  0 1 0 0 0 0 0 0;
                  0 0 1 0 0 0 0 0;
                  0 0 0 1 0 0 0 0;
                  0 0 0 0 0 0 0 1;
                  0 0 0 0 0 0 1 0;
                  0 0 0 0 0 1 0 0;
                  0 0 0 0 1 0 0 0];
        return;
    end
    return;
end
