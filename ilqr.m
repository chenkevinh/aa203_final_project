function [s_bar, u_bar, Y, y] = ilqr(As, Bs, s0, s_goal, N, Q, R, QN, eps, max_iters)

%     """Compute the iLQR set-point tracking solution.
%
%     Arguments
%     ---------
%     f : callable
%         A function describing the discrete-time dynamics, such that
%         `s(k+1) = f(s(k), u(k))`.
%     s0 : numpy.ndarray
%         The initial state (1-D).
%     s_goal : numpy.ndarray
%         The goal state (1-D).
%     N : int
%         The time horizon of the LQR cost function.
%     Q : numpy.ndarray
%         The state cost matrix (2-D).
%     R : numpy.ndarray
%         The control cost matrix (2-D).
%     QN : numpy.ndarray
%         The terminal state cost matrix (2-D).
%     eps : float, optional
%         Termination threshold for iLQR.
%     max_iters : int, optional
%         Maximum number of iLQR iterations.
%
%     Returns
%     -------
%     s_bar : numpy.ndarray
%         A 2-D array where `s_bar(k)` is the nominal state at time step `k`,
%         for `k = 0, 1, ..., N-1`
%     u_bar : numpy.ndarray
%         A 2-D array where `u_bar(k)` is the nominal control at time step `k`,
%         for `k = 0, 1, ..., N-1`
%     Y : numpy.ndarray
%         A 3-D array where `Y(:,:,k)` is the matrix gain term of the iLQR control
%         law at time step `k`, for `k = 0, 1, ..., N-1`
%     y : numpy.ndarray
%         A 2-D array where `y(:,k)` is the offset term of the iLQR control law
%         at time step `k`, for `k = 0, 1, ..., N-1`
%     """

n = size(Q,1);        % state dimension
m = size(R,1);        % control dimension

% Initialize gains `Y` and offsets `y` for the policy
Y = zeros(m, n, N);
y = zeros(m, N);

% Initialize the nominal trajectory `(s_bar, u_bar`), and the
% deviations `(ds, du)`
u_bar = zeros(m, N);
s_bar = zeros(n, N + 1);
s_bar(:,1) = s0;
for k = 1:N
    s_bar(:,k+1) = As(:,:,k) * s_bar(:,k) + Bs(:,:,k) * u_bar(:,k);
end
ds = zeros(n, N+1);
du = zeros(m, N);

% iLQR loop
converged = 0;
i = 0;
while i <= max_iters

    P = zeros(n,n,N+1);
    p = zeros(n,N+1);
    P(:,:,end) = QN;
    c_N = s_bar(:,end) - s_goal;
    p(:,end) = (0.5*c_N'*QN)';

    for k = N:-1:1
        c_k = s_bar(:,k) - s_goal;
        q_k = (0.5*c_k' * Q)';
        r_k = (0.5*u_bar(:,k)' * R)';
        h_xk = q_k + As(:,:,k)' * p(:,k+1);
        h_uk = r_k + Bs(:,:,k)' * p(:,k+1);
        Hxxk = Q + As(:,:,k)' * P(:,:,k+1) * As(:,:,k);
        Huuk = R + Bs(:,:,k)' * P(:,:,k+1) * Bs(:,:,k);
        Hxuk = As(:,:,k)' * P(:,:,k+1) * Bs(:,:,k);

        Y(:,:,k) = -Huuk^-1 * Hxuk';
        y(:,k) = -Huuk^-1 * h_uk;

        P(:,:,k) = Hxxk + Hxuk * Y(:,:,k);
        p(:,k) = h_xk + Hxuk * y(:,k);
    end

    for k = 1:N
        du(:,k) = y(:,k) + Y(:,:,k) * ds(:,k);
        ds(:,k+1) = As(:,:,k)*(s_bar(:,k) + ds(:,k)) + Bs(:,:,k) * (u_bar(:,k)+du(:,k)) - s_bar(:,k+1);
        s_bar(:,k) = s_bar(:,k) + ds(:,k);
        u_bar(:,k) = u_bar(:,k) + du(:,k);
    end

    if max(abs(du)) < eps
        converged = 1;
        break
    end
    i = i + 1;
end
if converged == 0
    disp('iLQR did not converge!')
end
end

