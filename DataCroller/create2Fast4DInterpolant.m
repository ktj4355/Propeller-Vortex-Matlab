
function F = create2Fast4DInterpolant(X,V)
% 한번만 수행하여 저장할 정보

ind_0=find(abs(X(:,1)-0)<0.01);
ind_25=find(abs(X(:,1)-0.25)<0.01);
ind_50=find(abs(X(:,1)-0.5)<0.01);
ind_75=find(abs(X(:,1)-0.75)<0.01);
ind_100=find((abs(X(:,1)-1)<0.01));

BLD0=scatteredInterpolant(X(ind_0,2),X(ind_0,3),X(ind_0,4),V(ind_0,1),"linear","nearest");
BLD25=scatteredInterpolant(X(ind_25,2),X(ind_25,3),X(ind_25,4),V(ind_25,1),"linear","nearest");
BLD50=scatteredInterpolant(X(ind_50,2),X(ind_50,3),X(ind_50,4),V(ind_50,1),"linear","nearest");
BLD75=scatteredInterpolant(X(ind_75,2),X(ind_75,3),X(ind_75,4),V(ind_75,1),"linear","nearest");
BLD100=scatteredInterpolant(X(ind_100,2),X(ind_100,3),X(ind_100,4),V(ind_100,1),"linear","nearest");

F=@(blend, thick, RE, AOA) interpFast(BLD0,BLD25,BLD50,BLD75,BLD100,[blend, thick, RE, AOA]);

    function Vq = interpFast(BLD0,BLD25,BLD50,BLD75,BLD100,queryPts)
        %queryPts [blend, thick, RE, AOA]

        D=floor(queryPts(1)/0.25);

        switch D
            case 0
                A=BLD0(queryPts(2),queryPts(3),queryPts(4));
                B=BLD25(queryPts(2),queryPts(3),queryPts(4));

                Vq=A+(B-A).*((queryPts(1)-0)./(0.25-0));
                return

            case 1
                A=BLD25(queryPts(2),queryPts(3),queryPts(4));
                B=BLD50(queryPts(2),queryPts(3),queryPts(4));

                Vq=A+(B-A).*((queryPts(1)-0.25)./(0.5-0.25));
                return
            case 2
                A=BLD50(queryPts(2),queryPts(3),queryPts(4));
                B=BLD75(queryPts(2),queryPts(3),queryPts(4));

                Vq=A+(B-A).*((queryPts(1)-0.5)./(0.75-0.5));
                return
            case 3
                A=BLD75(queryPts(2),queryPts(3),queryPts(4));
                B=BLD100(queryPts(2),queryPts(3),queryPts(4));

                Vq=A+(B-A).*((queryPts(1)-0.75)./(1-0.75));
                return
            case 4

                B=BLD100(queryPts(2),queryPts(3),queryPts(4));

                Vq=B;
                return

        end

    end


end

