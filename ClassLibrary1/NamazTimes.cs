using System;
using System.Collections.Generic;
using System.Text;
namespace Islam
{

    public enum CalculationMehod
    {
        Jafari, Karachi, ISNA, MWL, Makkah, Egypt, Custom, Tehran
    }
    public enum JuristicMethod
    {
        Shafii = 1, Hanafi
    }
    public enum HigherLattitude
    {
        None, MidNight, OneSeventh, AngleBased
    }
    public enum TimeFormat
    {
        Time24, Time12
    }

    public class NamazTimes
    {
        private double lattitude;
        private double longitude;
        private double timeZone;
        private DateTime dateTime;

        private double JDate;
        private double[][] methodParams;


        public NamazTimes()
        {
            this.calcMethod = CalculationMehod.Karachi;
            this.highLats = HigherLattitude.None;
            this.timeFormat = TimeFormat.Time12;
            this.juristicMethod = JuristicMethod.Hanafi;

            methodParams = new double[][]
            {
                new double[]{16, 0, 4, 0, 14},
                new double[]{18, 1, 0, 0, 18},
                new double[]{15, 1, 0, 0, 15},
                new double[]{18, 1, 0, 0, 17},
                new double[]{19, 1, 0, 1, 90},
                new double[]{19.5, 1, 0, 0, 17.5},
                new double[]{17.7, 0, 4.5, 0, 15},
                new double[]{18, 1, 0, 0, 17}
            };
        }

        public DateTime[] GetPrayerTimes()
        {
            return GetPrayerTimes(this.dateTime, this.lattitude, this.longitude, this.timeZone);
        }
        public int[] GetPrayerMinutes()
        {
            int[] times = new int[8];
            DateTime[] dt = GetPrayerTimes();
            TimeSpan ts = new TimeSpan();
            for (int i = 0; i <= 5; i++)
            {
                ts = dt[i + 1].TimeOfDay - dt[i].TimeOfDay;
                times[i] = (int)ts.TotalMinutes;
            }

            //ts = dt[2].TimeOfDay - dt[1].TimeOfDay;
            //times[1] = (int)ts.TotalMinutes;
            //ts = dt[3].TimeOfDay - dt[2].TimeOfDay;
            //times[2] = (int)ts.TotalMinutes;
            //ts = dt[4].TimeOfDay - dt[3].TimeOfDay;
            //times[3] = (int)ts.TotalMinutes;
            //ts = dt[5].TimeOfDay - dt[4].TimeOfDay;
            //times[4] = (int)ts.TotalMinutes;
            //ts = dt[6].TimeOfDay - dt[5].TimeOfDay;
            //times[5] = (int)ts.TotalMinutes;

            //isha to 12 am minutes are calculated below
            ts = new DateTime(DateTime.Now.Year, DateTime.Now.Month, DateTime.Now.Day, 23, 59, 59) - DateTime.Parse(GetPrayerTimes()[6].ToShortTimeString());
            times[6] = (int)ts.TotalMinutes;//TotalNightMinutes - (times[0] + times[5]);
            times[7] = (int)GetPrayerTimes()[0].TimeOfDay.TotalMinutes;
            return times;
        }
        public DateTime[] GetPrayerTimes(DateTime date, double lattitude, double longitude, double timeZone)
        {
            this.lattitude = lattitude;
            this.longitude = longitude;
            this.timeZone = timeZone;
            this.JDate = this.JulianDate(date) - longitude / (15 * 24);
            double[] times = this.ComputeTimes();
            DateTime[] dt = new DateTime[7];
            times = this.AdjustTimes(times);
            string[] tempStr = this.AdjustTimesFormat(times);
            //SalaahTimeCalculator stc = new SalaahTimeCalculator();
            //string[] tempStr = stc.getPrayerTimes(date, lattitude, longitude, (int)timeZone);
            for (int i = 0; i <= 6; i++)
            {
                dt[i] = DateTime.Parse(tempStr[i]);
            }
            return dt;
        }

        public DateTime GetTime(int TimeNumber)
        {
            //string outStr="";
            DateTime[] Times = GetPrayerTimes(this.dateTime, this.lattitude, this.longitude, this.timeZone);
            switch (TimeNumber)
            {
                case 0:
                    return Times[0];
                case 1:
                    return Times[1];

                case 2:
                    return Times[2];

                case 3:
                    return Times[3];

                case 4:
                    return Times[4];

                case 5:
                    return Times[5];

                case 6:
                    return Times[6];
                default:
                    return Times[0];

            }

            //return DateTime.Parse(outStr);
        }
        public double JulianDate(DateTime date)
        {
            int year = date.Year, month = date.Month, day = date.Day;

            if (date.Month <= 2)
            {
                year -= 1;
                month += 12;
            }

            double a = Math.Floor((double)year / 100);
            double b = 2 - a + Math.Floor(a / 4);
            return Math.Floor(365.25 * (year + 4716)) + Math.Floor(30.6001 * (month + 1)) + day + b - 1524.5;
        }

        double[] ComputeTimes()
        {
            double[] t = new double[] { 5, 6, 12, 13, 18, 18, 18 };
            for (int i = 0; i < 7; i++)
                t[i] /= 24;
            double Fajr = this.ComputeTime(180 - this.methodParams[(int)calcMethod][0], t[0]);
            double Sunrise = this.ComputeTime(180 - 0.833, t[1]);
            double Dhuhr = this.ComputeMidDay(t[2]);
            double Asr = this.ComputeAsr(this.juristicMethod, t[3]);
            double Sunset = this.ComputeTime(0.833, t[4]); ;
            double Maghrib = this.ComputeTime(this.methodParams[(int)calcMethod][2], t[5]);
            double Isha = this.ComputeTime(this.methodParams[(int)calcMethod][4], t[6]);

            return new double[] { Fajr, Sunrise, Dhuhr, Asr, Sunset, Maghrib, Isha };
        }

        double ComputeTime(double angle, double dayPortion)
        {
            double d = this.SunDeclination(this.JDate + dayPortion);
            double z = this.ComputeMidDay(dayPortion);
            double a = (-this.Dsin(angle) - this.Dsin(d) * this.Dsin(this.lattitude)) /
                    (this.Dcos(d) * this.Dcos(this.lattitude));
            double v = 1d / 15 * this.Darccos(a);
            return z + (angle > 90 ? -v : v);
        }

        public double SunDeclination(double julianDate)
        {
            double d = julianDate - 2451545.0;
            double g = this.Fixangle(357.529 + 0.98560028 * d);
            double q = this.Fixangle(280.459 + 0.98564736 * d);
            double l = this.Fixangle(q + 1.915 * this.Dsin(g) + 0.020 * this.Dsin(2 * g));
            double e = 23.439 - 0.00000036 * d;
            return this.Darcsin(this.Dsin(e) * this.Dsin(l));
        }

        public double EquationOfTime(double julianDate)
        {
            double d = julianDate - 2451545.0;
            double g = this.Fixangle(357.529 + 0.98560028 * d);
            double q = this.Fixangle(280.459 + 0.98564736 * d);
            double l = this.Fixangle(q + 1.915 * this.Dsin(g) + 0.020 * this.Dsin(2 * g));
            double e = 23.439 - 0.00000036 * d;
            double ra = this.Darctan2(this.Dcos(e) * this.Dsin(l), this.Dcos(l)) / 15;
            ra = this.Fixhour(ra);
            return q / 15 - ra;
        }

        double ComputeMidDay(double dayPortion)
        {
            double t = this.EquationOfTime(this.JDate + dayPortion);
            return this.Fixhour(12 - t);
        }
        double ComputeAsr(JuristicMethod juristicMethod, double dayPortion)  // Shafii: step=1, Hanafi: step=2
        {
            double d = this.SunDeclination(this.JDate + dayPortion);
            double g = -this.Darccot(((double)juristicMethod) + this.Dtan(Math.Abs(this.lattitude - d)));
            return this.ComputeTime(g, dayPortion);
        }

        double[] AdjustTimes(double[] times)
        {
            for (int i = 0; i < 7; i++)
                times[i] += this.timeZone - this.longitude / 15;
            times[2] += this.dhuhrMinutes / 60; //Dhuhr
            if (this.methodParams[(int)this.calcMethod][1] == 1) // Maghrib
                times[5] = times[4] + this.methodParams[(int)this.calcMethod][2] / 60;
            if (this.methodParams[(int)this.calcMethod][3] == 1) // Isha
                times[6] = times[5] + this.methodParams[(int)this.calcMethod][4] / 60;

            if (this.highLats != HigherLattitude.None)
                times = this.AdjustHighLatTimes(times);
            return times;
        }

        string[] AdjustTimesFormat(double[] times)
        {
            string[] temp = new string[times.Length];

            for (int i = 0; i < 7; i++)
            {
                if (this.timeFormat == TimeFormat.Time12)
                    temp[i] = this.FloatToTime12(times[i], false);
                else if (this.timeFormat == TimeFormat.Time24)
                    temp[i] = this.FloatToTime24(times[i]);
            }
            return temp;
        }

        double[] AdjustHighLatTimes(double[] times)
        {
            double nightTime = this.TimeDiff(times[4], times[1]); // sunset to sunrise

            // Adjust Fajr
            double FajrDiff = this.NightPortion(this.methodParams[(int)this.calcMethod][0]) * nightTime;
            if (times[0] == 0d || this.TimeDiff(times[0], times[1]) > FajrDiff)
                times[0] = times[1] - FajrDiff;

            // Adjust Isha
            double IshaAngle = (this.methodParams[(int)this.calcMethod][3] == 0) ? this.methodParams[(int)this.calcMethod][4] : 18;
            double IshaDiff = this.NightPortion(IshaAngle) * nightTime;
            if (times[6] == 0d || this.TimeDiff(times[4], times[6]) > IshaDiff)
                times[6] = times[4] + IshaDiff;

            // Adjust Maghrib
            double MaghribAngle = (this.methodParams[(int)this.calcMethod][1] == 0) ? this.methodParams[(int)this.calcMethod][2] : 4;
            double MaghribDiff = this.NightPortion(MaghribAngle) * nightTime;
            if (times[5] == 0d || this.TimeDiff(times[4], times[5]) > MaghribDiff)
                times[5] = times[4] + MaghribDiff;

            return times;
        }

        double NightPortion(double angle)
        {
            double r = 0d;
            if (this.highLats == HigherLattitude.AngleBased)
                r = 1d / 60 * angle;
            if (this.highLats == HigherLattitude.MidNight)
                r = 1d / 2;
            if (this.highLats == HigherLattitude.OneSeventh)
                r = 1d / 7;
            return r;
        }

        double TimeDiff(double time1, double time2)
        {
            return this.Fixhour(time2 - time1);
        }


        string FloatToTime12(double time, bool noSuffix)
        {
            if (time == 0d)
                throw new Exception("Invalid Time");
            time = this.Fixhour(time + 0.5 / 60);  // add 0.5 minutes to round
            double hours = Math.Floor(time);
            double minutes = Math.Floor((time - hours) * 60);
            string suffix = hours >= 12 ? " pm" : " am";
            hours = (hours + 12 - 1) % 12 + 1;
            return hours + ":" + this.TwoDigitsFormat(minutes) + (noSuffix ? "" : suffix.ToString());
        }

        string FloatToTime24(double time)
        {
            if (time == 0d)
                throw new Exception("Invalid Time");
            time = this.Fixhour(time + 0.5 / 60);  // add 0.5 minutes to round
            double hours = Math.Floor(time);
            double minutes = Math.Floor((time - hours) * 60);
            return this.TwoDigitsFormat(hours) + ":" + this.TwoDigitsFormat(minutes);
        }

        string TwoDigitsFormat(double num)
        {
            return (num < 10) ? "0" + num : num.ToString(); ;
        }

        double Darctan2(double y, double x)
        {
            return this.Dtd(Math.Atan2(y, x));
        }

        double Fixhour(double a)
        {
            a = a - 24.0 * (Math.Floor(a / 24.0));
            return a < 0 ? a + 24.0 : a;
        }

        double Fixangle(double a)
        {
            a = a - 360.0 * (Math.Floor(a / 360.0));
            return a < 0 ? a + 360.0 : a;
        }
        double Darccos(double x)
        {
            return this.Dtd(Math.Acos(x));
        }
        double Darcsin(double x)
        {
            return this.Dtd(Math.Asin(x));
        }
        double Darccot(double x)
        {
            return this.Dtd(Math.Atan(1d / x));
        }

        double Dsin(double d)
        {
            return Math.Sin(this.Dtr(d));
        }
        double Dcos(double d)
        {
            return Math.Cos(this.Dtr(d));
        }
        double Dtan(double d)
        {
            return Math.Tan(this.Dtr(d));
        }
        double Dtd(double r)
        {
            return (r * 180.0) / Math.PI;
        }

        double Dtr(double r)
        {
            return (r * Math.PI) / 180.0;
        }

        int GetMinutes(DateTime time1, DateTime time2, bool complement)
        {
            TimeSpan ts = new TimeSpan();
            TimeSpan ts24 = new TimeSpan(24, 0, 0);
            if (complement)
                ts = ts24 + (time2 - time1);
            else
                ts = time2 - time1;
            return Math.Abs((int)ts.TotalMinutes);
        }
        int TDayMin(string sunrise, string sunset)
        {
            DateTime Sunrise = DateTime.Parse(sunrise);
            DateTime Sunset = DateTime.Parse(sunset);
            TimeSpan ts = new TimeSpan();
            ts = Sunset - Sunrise;
            return Math.Abs((int)ts.TotalMinutes);
        }
        int PDayMin(string sunset, string sunrise)
        {
            DateTime Sunset = DateTime.Parse(sunset);
            DateTime Sunrise = DateTime.Parse(sunrise);
            DateTime Nowtime = DateTime.Parse(DateTime.Now.ToShortTimeString());
            TimeSpan ts = new TimeSpan();
            if (Nowtime > Sunrise && Nowtime < Sunset)
            {
                ts = Nowtime.TimeOfDay - Sunrise.TimeOfDay;
                return Math.Abs((int)ts.TotalMinutes);
            }
            else
                return TDayMin(sunrise, sunset);


        }
        int TNightMin(string sunrise, string sunset)
        {

            DateTime Sunrise = DateTime.Parse(sunrise);
            DateTime Sunset = DateTime.Parse(sunset);
            TimeSpan ts = new TimeSpan();
            ts = new TimeSpan(24, 0, 0) - Sunset.TimeOfDay + Sunrise.TimeOfDay;
            return Math.Abs((int)ts.TotalMinutes);
        }
        int PNightMin(string sunset, string sunrise)
        {
            DateTime Nowtime = DateTime.Parse(DateTime.Now.ToShortTimeString());
            DateTime Sunset = DateTime.Parse(sunset);
            DateTime Sunrise = DateTime.Parse(sunrise);
            TimeSpan ts = new TimeSpan();
            DateTime zero00 = DateTime.Parse("23:59:59 PM");

            if (Nowtime < zero00 && Nowtime > Sunset)
            {
                ts = Nowtime - Sunset;
                return Math.Abs((int)ts.TotalMinutes);
            }
            else if (Nowtime < Sunrise)
            {
                ts = zero00 - Sunset;
                ts += Nowtime.TimeOfDay;
                return (int)ts.TotalMinutes;
            }
            else
                return 0;
        }

        public int TotalDayMinutes
        {
            get { return TDayMin(GetTime(1).ToShortTimeString(), GetTime(4).ToShortTimeString()); }
        }

        public int PassedDayMintues
        {
            get { return PDayMin(GetTime(4).ToShortTimeString(), GetTime(1).ToShortTimeString()); }
        }
        public int TotalNightMinutes
        {
            get { return TNightMin(GetTime(1).ToShortTimeString(), GetTime(4).ToShortTimeString()); }
        }
        public int PassedNightMinutes
        {
            get { return PNightMin(GetTime(4).ToShortTimeString(), GetTime(1).ToShortTimeString()); }
        }


        private TimeFormat timeFormat;

        public TimeFormat TimeFormat
        {
            get { return timeFormat; }
            set { timeFormat = value; }
        }

        public double Lattitude
        {
            get { return this.lattitude; }
            set { this.lattitude = value; }
        }
        public double Longitude
        {
            get { return this.longitude; }
            set { this.longitude = value; }
        }
        public double Zone
        {
            get { return this.timeZone; }
            set { this.timeZone = value; }
        }

        public DateTime Date
        {
            get { return this.dateTime; }
            set { this.dateTime = value; }
        }





        private CalculationMehod calcMethod;

        public CalculationMehod CalcMethod
        {
            get { return calcMethod; }
            set { calcMethod = value; }
        }

        private JuristicMethod juristicMethod;

        public JuristicMethod JuristicMethod
        {
            get { return juristicMethod; }
            set { juristicMethod = value; }
        }

        private double dhuhrMinutes;

        public double DhuhrMinutes
        {
            get { return dhuhrMinutes; }
            set { dhuhrMinutes = value; }
        }

        private HigherLattitude highLats;

        public HigherLattitude HigherLattitude
        {
            get { return highLats; }
            set { highLats = value; }
        }

        public string ToString(DateTime date, double lattitude, double longitude, double timeZone)
        {
            DateTime[] strs = GetPrayerTimes(date, lattitude, longitude, timeZone);
            string result = "Date : " + date.ToShortDateString() + "\n";

            result += "Fajr : " + strs[0].ToShortTimeString() + "\n";
            result += "Sunrise : " + strs[1].ToShortTimeString() + "\n";
            result += "Dhuhr : " + strs[2].ToShortTimeString() + "\n";
            result += "Asr : " + strs[3].ToShortTimeString() + "\n";
            result += "Sunset : " + strs[4].ToShortTimeString() + "\n";
            result += "Maghrib : " + strs[5].ToShortTimeString() + "\n";
            result += "Isha : " + strs[6].ToShortTimeString() + "\n";

            return result;
        }
    }
}

