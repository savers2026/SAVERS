% Amir

function ca = generate_CA_code(prn)
    % Generate one GPS C/A PRN
    g1 = -ones(1, 1023);
    g2 = -ones(1, 1023);

    % G2 tap mapping for all PRNs
    tap = [2 6;3 7;4 8;5 9;1 9;2 10;1 8;2 9;3 10;2 3; ...
           3 4;5 6;6 7;7 8;8 9;9 10;1 4;2 5;3 6;4 7; ...
           5 8;6 9;1 3;4 6;5 7;6 8;7 9;8 10;1 6;2 7; ...
           3 8;4 9;5 10;4 10;1 7;2 8;4 10;1 8;2 9;3 10;4 10];

    t1 = tap(prn,1);
    t2 = tap(prn,2);

    % shift registers
    R1 = ones(1,10);
    R2 = ones(1,10);

    for i = 1:1023
        g1(i) = R1(10);
        g2(i) = xor(R2(t1), R2(t2));

        R1 = [xor(R1(3), R1(10)), R1(1:9)];
        temp = xor(R2(2), xor(R2(3), xor(R2(6), xor(R2(8), R2(9)))));
        R2 = [temp, R2(1:9)];
    end

    % C/A = G1 ⊕ G2 → convert to ±1
    ca = 1 - 2 * xor(g1, g2);
end
