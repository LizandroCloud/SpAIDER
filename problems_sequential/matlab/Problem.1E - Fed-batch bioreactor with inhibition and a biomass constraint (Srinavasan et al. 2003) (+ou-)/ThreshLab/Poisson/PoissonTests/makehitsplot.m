if exist('hits') ~= 1,
   hits20031205
end
if ~exist('plotoption'),
   plotoption = '-';
end
N = length(hits);
plot(hits,plotoption);
week1 = 5;
weeks = (week1:week1+N-1);
week1 = 1;
year = 1997;
firstdayofyear = 3; % Sunday = 0;
while week1 < N,
   daysinyear = 365 + (mod(year,4) == 0) - (mod(year,100) == 0) ...
                                         + (mod(year,400) == 0);
   lastdayofyear = firstdayofyear + mod(daysinyear,7) - 1;
   lastdayofyear = mod(lastdayofyear,7);
   daysinyear = daysinyear + 6 - lastdayofyear;
   daysinyear = daysinyear - mod(7 - firstdayofyear,7);
   % daysinyear is now number of days from first sunday of year until first
   % sunday of next year.
   weeksinyear = div(daysinyear,7);
   weeks(week1:N) = weeks(week1:N) - weeksinyear*(weeks(week1:N)>weeksinyear);

   year = year + 1;
   week1 = week1+weeksinyear;
   firstdayofyear = mod(lastdayofyear + 1,7);
end

lastweeks = zeros(size(weeks));
lastweeks(1:N-1) = (weeks(2:N) == 1);
hmax = max(hits);
hold on
lastweeks = find(lastweeks);
for week=lastweeks,
   plot([week week],[0,hmax],'b:')
end
hold off
axis([1,N,0,hmax])

lastweeks = [0 lastweeks N];
r = 0.7;
for y = 1:year-1-1996,
   text(r*lastweeks(y)+(1-r)*lastweeks(y+1),0.9*hmax,num2str(1996+y));
end

text(N-30,0.025*hmax,'weeks')
text(5,0.975*hmax,'# hits')
[sorthits i] = sort(hits);
rang = (N:-1:1);
rang(i) = rang;
% !! This is not the same as rang = rang(i) !!
yp0 = hits;
% text(12,1400,'1')
% text(170,6500,'2')
