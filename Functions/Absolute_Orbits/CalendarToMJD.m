function MJD = CalendarToMJD(UT1)
    % CALENDARTOMJD takes in UT1 time in calendar format and converts it to
    % MJD format. 
    % The input UT1 time must be in format [YYYY, MM, DD].
    % Note that DD is represented as an integer representing the day of the
    % month plus the hours that have passed in that day represented as a
    % fractional number of days.
    % i.e. Jan 1st 12:00 noon would have DD = 1.5 even though technically
    % only 0.5 days have occured during that month at that point.
    Y = UT1(1);
    M = UT1(2);
    D = UT1(3);
    y = Y;
    m = M;
    if M <= 2
        y = y - 1;
        m = m + 12;
    end
    % B = floor(Y/400) - floor(Y/100) + floor(Y/4);
    B = Y/400 - Y/100 + Y/4;
    MJD = (365*y) - 679004 + floor(B) + floor(30.6001*(m+1)) + D;
end
