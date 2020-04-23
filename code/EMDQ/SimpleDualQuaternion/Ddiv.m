function [ dn ] = Ddiv( dn0, dn1 )
    dn = zeros(size(dn0,1),2);
	dn(:,1) = dn0(:,1) ./ dn1(:,1);
	dn(:,2) = dn0(:,2) ./ dn1(:,1) - dn0(:,1) * dn1(:,2) ./ (dn1(:,1) * dn1(:,1));

end

