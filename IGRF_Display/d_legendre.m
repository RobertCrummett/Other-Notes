%% Source material
% On the Computation af Derivatives of Legendre Functions, W. Bosch, 2000
% DOI:10.1016/S1464-1895(00)00101-0

%% P = D_LEGENDRE(N,X,NORMALIZE)
% Function behaves exactly like MATLAB's LEGENDRE function
% See LEGENDRE for details on inputs

function P = d_legendre(N, X, normalize)
% Better than the online forum function because it is stable at poles
% Last edit 4/2/2023 Nate

% Compute P_{nm}
C = legendre(N, X, normalize);
[MLEN,XLEN] = size(C);
P = zeros(MLEN,XLEN);

% Compute d(P_{nm})/d(theta)
switch normalize

    % Unnormalized functions
    case 'unnorm'
        
        % Order zero
        if N == 0
            return % (10), d(P_{00})/d(theta) := 0
        
        % Order one or more
        else
            P(1,:) = -C(2,:); % (9)
            MIDX = 2:(MLEN-1); M = repmat(MIDX',1,XLEN) - 1;
            P(2:(MLEN-1),:) = (N + M).*(N - M + 1).*C(MIDX-1,:)/2 - ...
                C(MIDX+1,:)/2; % (8)
            P(MLEN,:) = N*C(MLEN-1,:); % (11)
            return
        end

    % Normalized functions
    otherwise

        % Order zero
        if N == 0
            return % A7, d(P_{00})/d(theta) := 0

        % Order one
        elseif N == 1
            P(1,:) = -C(2,:); % A7, N = 1
            P(2,:) = C(1,:); % A8
            return

        % Order two or more
        else
            P(1,:) = -sqrt(N*(N+1)/2)*C(2,:); % A7
            P(2,:) = sqrt(N*(N+1))*C(1,:)/2 - sqrt((N-1)*(N+2))*C(3,:)/2; % A8

            % Order two
            if N == 2
                P(3,:) = C(2,:); % A9
                return
            else
                MIDX = 3:(MLEN-1); M = repmat(MIDX',1,XLEN) - 1;
                P(MIDX,:) = sqrt((N + M).*(N - M + 1)).*C(MIDX-1,:)/2-...
                    sqrt((N + M + 1).*(N - M)).*C(MIDX+1,:)/2; % A6
                P(MLEN,:) = sqrt(N/2)*C(MLEN-1,:); % A9
                return
            end
        end
end

