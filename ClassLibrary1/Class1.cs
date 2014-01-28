using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;

public enum CalculationMethods
{
    Jafari,
    Karachi,
    ISNA,
    MWL,
    Makkah,
    Egypt,
    Custom
}

public enum JurusticMethods
{
    Shafii,
    Hanafi
}

namespace Islam
{
    /// <summary>
    /// Current Author : Jameel Haffejee
    /// Contact Info   : JameelHaffejee@gmail.com
    /// Description    :
    /// A salaah time calculator. This class calculates the time based on the users lattitude and longitude and their current timezone.
    /// Note that this is a direct conversion from the javascript version.
    /// Please note that even though i did not design the original code i did do the
    /// entire conversion my self , so please credit where needed.
    /// Also if possible let me know that you are using the code and wether you are happy with it or not.
    /// Modifications had to be made as needed to allow the conversion to .Net .
    /// Original  version here : http://tanzil.info/praytime/doc/manual/
    ///
    /// The License below is included as per the original code.
    /// 
    ///PrayTime: Prayer Times Calculator (ver 1.1)
    ///Copyright (C) 2007, Hamid Zarrabi-Zadeh
    ///This program is free software; you can redistribute it and/or
    ///modify it under the terms of the GNU General Public License
    ///as published by the Free Software Foundation; either version 2
    ///of the License, or (at your option) any later version.
    ///This program is distributed in the hope that it will be useful,
    ///but WITHOUT ANY WARRANTY; without even the implied warranty of
    ///MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    ///GNU General Public License for more details.
    ///You can get a copy of the GNU General Public License from
    ///http://www.gnu.org/copyleft/gpl.html
    /// </summary>
    public class SalaahTimeCalculator
    {
        #region Constants
        // Adjusting Methods for Higher Latitudes
        int None = 0;    // No adjustment
        int MidNight = 1;    // middle of night
        int OneSeventh = 2;    // 1/7th of night
        int AngleBased = 3;    // angle/60th of night


        // Time Formats
        int Time24 = 0;    // 24-hour format
        int Time12 = 1;    // 12-hour format
        int Time12NS = 2;    // 12-hour format with no suffix
        int Float = 3;    // floating point number 

        // Time Names
        string[] timeNames = new string[]{
		"Fajr",
		"Sunrise",
		"Dhuhr",
		"Asr",
		"Sunset",
		"Maghrib",
		"Isha"
	};

        #endregion

        #region Global doubleiables
        int calcMethod = 0;		// caculation method
        int asrJuristic = 0;		// Juristic method for Asr
        int dhuhrMinutes = 0;		// minutes after mid-day for Dhuhr
        int adjustHighLats = 1;	// adjusting method for higher latitudes

        int timeFormat = 0;		// time format

        double lat;        // latitude 
        double lng;        // longitude 
        int timeZone;   // time-zone 
        double JDate;      // Julian date
        #endregion

        #region Technical Settings
        int numIterations = 1;		// number of iterations needed to compute times, this should never be more than 1;
        #endregion

        #region Calc Method Parameters
        /// <summary>
        ///  methodParams[methodNum] = new Array(fa, ms, mv, is, iv);	
        ///     fa : fajr angle
        ///     ms : maghrib selector (0 = angle; 1 = minutes after sunset)
        ///     mv : maghrib parameter value (in angle or minutes)
        ///     is : isha selector (0 = angle; 1 = minutes after maghrib)
        ///     iv : isha parameter value (in angle or minutes)
        /// </summary>
        double[][] methodParams = new double[][]{
                    new double[]{16, 0, 4, 0, 14},//Jafari
                    new double[]{18, 1, 0, 0, 18},//Karachi
                    new double[]{15, 1, 0, 0, 15},//ISNA
                    new double[]{18, 1, 0, 0, 17},//MWL
                    new double[]{19, 1, 0, 1, 90},//Makkah
                    new double[]{19.5, 1, 0, 0, 17.5},//Egypt
                    new double[]{18, 1, 0, 0, 17}//Custom
                        };

        #endregion

        #region Interface Functions

        ///<summary>
        /// Returns the prayer times for a given date , the date format is specified as indiviual settings.
        /// </summary>
        /// <param name="year">Year to use when calculating times</param>        
        /// <param name="month">Month to use when calculating times</param>
        /// <param name="day">Day to use when calculating times</param>
        /// <param name="latitude">Lattitude to use when calculating times</param>
        /// <param name="longitude">Longitude to use when calculating times</param>
        /// <param name="timeZone">Time zone to use when calculating times</param>
        /// <returns>
        /// A string Array containing the Salaah times
        /// </returns>
        private string[] getDatePrayerTimes(int year, int month, int day, double latitude, double longitude, int timeZone)
        {
            lat = latitude;
            lng = longitude;
            this.timeZone = effectiveTimeZone(year, month, day, timeZone);
            JDate = julianDate(year, month, day) - longitude / (15 * 24);
            return computeDayTimes();
        }

        ///<summary>
        /// Returns the prayer times for a given date , the date format is specified as indiviual settings.
        /// </summary>
        /// <param name="year">Year to use when calculating times</param>        
        /// <param name="month">Month to use when calculating times</param>
        /// <param name="day">Day to use when calculating times</param>
        /// <param name="latitude">Lattitude to use when calculating times</param>
        /// <param name="longitude">Longitude to use when calculating times</param>
        /// <param name="timeZone">Time zone to use when calculating times</param>
        /// <returns>
        /// A string Array containing the Salaah times,The time is in the 24 hour format.
        /// The array is structured as such.
        /// 0.Fajr
        /// 1.Sunrise
        /// 2.Dhuhr
        /// 3.Asr
        /// 4.Sunset
        /// 5.Maghrib
        /// 6.Isha
        /// </returns>
        public string[] getPrayerTimes(DateTime date, double latitude, double longitude, int timeZone)
        {
            return this.getDatePrayerTimes(date.Year, date.Month, date.Day,
                        latitude, longitude, timeZone);
        }


        ///<summary>
        ///This method is used to set the calculation method for the salaah times.
        ///There are 7 calculation types available.
        ///Jafari  = 0
        ///Karachi = 1 
        ///ISNA    = 2
        ///MWL     = 3
        ///Makkah  = 4
        ///Egypt   = 5
        ///Custom  = 6
        /// </summary>
        /// <param name="methodToUse">Calculation Method to use</param>
        public void setCalculationMethod(CalculationMethods methodToUse)
        {
            this.calcMethod = (int)methodToUse;
        }

        // set the juristic method for Asr
        public void setAsrJurusticionType(JurusticMethods selectedJuristicion)
        {
            this.asrJuristic = (int)selectedJuristicion;
        }

        // set the angle for calculating Fajr
        private void setFajrAngle(double angle)
        {
            this.setCustomParams(new double?[] { angle, null, null, null, null });
        }

        // set the angle for calculating Maghrib
        private void setMaghribAngle(double angle)
        {
            this.setCustomParams(new double?[] { null, 0, angle, null, null });
        }

        // set the angle for calculating Isha
        private void setIshaAngle(double angle)
        {
            this.setCustomParams(new double?[] { null, null, null, 0, angle });
        }


        // set the minutes after mid-day for calculating Dhuhr
        private void setDhuhrMinutes(int minutes)
        {
            this.dhuhrMinutes = minutes;
        }

        // set the minutes after Sunset for calculating Maghrib
        private void setMaghribMinutes(int minutes)
        {
            this.setCustomParams(new double?[] { null, 1, minutes, null, null });
        }

        // set the minutes after Maghrib for calculating Isha
        private void setIshaMinutes(int minutes)
        {
            this.setCustomParams(new double?[] { null, null, null, 1, minutes });
        }

        // set custom values for calculation parameters
        private void setCustomParams(double?[] userParams)
        {
            for (int i = 0; i < 5; i++)
            {
                if (userParams[i] == null)
                    this.methodParams[(int)CalculationMethods.Custom][i] = this.methodParams[this.calcMethod][i];
                else
                    this.methodParams[(int)CalculationMethods.Custom][i] = userParams[i].Value;
            }
            this.calcMethod = (int)CalculationMethods.Custom;
        }

        // set adjusting method for higher latitudes 
        private void setHighLatsMethod(int methodID)
        {
            this.adjustHighLats = methodID;
        }

        // set the time format 
        public void setTimeFormat(int timeFormat)
        {
            this.timeFormat = timeFormat;
        }

        // convert float hours to 24h format
        private string floatToTime24(double time)
        {
            if (double.IsNaN(time))
                return "";
            time = this.fixhour(time + 0.5 / 60);  // add 0.5 minutes to round
            double hours = Math.Floor(time);
            double minutes = Math.Floor((time - hours) * 60);
            return this.twoDigitsFormat(Convert.ToInt32(hours)) + ":" + this.twoDigitsFormat(Convert.ToInt32(minutes));
        }

        // convert float hours to 12h format
        private string floatToTime12(double time, bool noSuffix)
        {
            if (double.IsNaN(time))
                return "";
            time = this.fixhour(time + 0.5 / 60);  // add 0.5 minutes to round
            double hours = Math.Floor(time);
            double minutes = Math.Floor((time - hours) * 60);
            string suffix = hours >= 12 ? " pm" : " am";
            hours = (hours + 12 - 1) % 12 + 1;
            return hours + ":" + this.twoDigitsFormat(Convert.ToInt32(minutes)) + (noSuffix ? "" : suffix);
        }

        // convert float hours to 12h format
        private string floatToTime12(double time)
        {
            bool noSuffix = false;
            if (double.IsNaN(time))
                return "";
            time = this.fixhour(time + 0.5 / 60);  // add 0.5 minutes to round
            double hours = Math.Floor(time);
            double minutes = Math.Floor((time - hours) * 60);
            string suffix = hours >= 12 ? " pm" : " am";
            hours = (hours + 12 - 1) % 12 + 1;
            return hours + ":" + this.twoDigitsFormat(Convert.ToInt32(minutes)) + (noSuffix ? "" : suffix);
        }

        // convert float hours to 12h format with no suffix
        private string floatToTime12NS(double time)
        {
            return this.floatToTime12(time, true);
        }
        #endregion

        #region Calculation Functions
        // References:
        // http://www.ummah.net/astronomy/saltime  
        // http://aa.usno.navy.mil/faq/docs/SunApprox.html


        // compute declination angle of sun and equation of time
        private double[] sunPosition(double jd)
        {
            double D = jd - 2451545.0;
            double g = this.fixangle(357.529 + 0.98560028 * D);
            double q = this.fixangle(280.459 + 0.98564736 * D);
            double L = this.fixangle(q + 1.915 * this.dsin(g) + 0.020 * this.dsin(2 * g));

            double R = 1.00014 - 0.01671 * this.dcos(g) - 0.00014 * this.dcos(2 * g);
            double e = 23.439 - 0.00000036 * D;

            double d = this.darcsin(this.dsin(e) * this.dsin(L));
            double RA = this.darctan2(this.dcos(e) * this.dsin(L), this.dcos(L)) / 15;
            RA = this.fixhour(RA);
            double EqT = q / 15 - RA;

            return new double[] { d, EqT };
        }

        // compute equation of time
        private double equationOfTime(double jd)
        {
            return this.sunPosition(jd)[1];
        }

        // compute declination angle of sun
        private double sunDeclination(double jd)
        {
            return this.sunPosition(jd)[0];
        }

        // compute mid-day (Dhuhr, Zawal) time
        private double computeMidDay(double t)
        {
            double T = this.equationOfTime(this.JDate + t);
            double Z = this.fixhour(12 - T);
            return Z;
        }

        // compute time for a given angle G
        private double computeTime(double G, double t)
        {
            double D = this.sunDeclination(this.JDate + t);
            double Z = this.computeMidDay(t);
            double V = ((double)(1 / 15d)) * this.darccos((-this.dsin(G) - this.dsin(D) * this.dsin(this.lat)) /
                    (this.dcos(D) * this.dcos(this.lat)));
            return Z + (G > 90 ? -V : V);
        }

        // compute the time of Asr
        private double computeAsr(int step, double t)  // Shafii: step=1, Hanafi: step=2
        {
            double D = this.sunDeclination(this.JDate + t);
            double G = -this.darccot(step + this.dtan(Math.Abs(this.lat - D)));
            return this.computeTime(G, t);
        }

        #endregion

        #region Compute Prayer Times
        // compute prayer times at given julian date
        private double[] computeTimes(double[] times)
        {
            double[] t = this.dayPortion(times);

            double Fajr = this.computeTime(180 - this.methodParams[this.calcMethod][0], t[0]);
            double Sunrise = this.computeTime(180 - 0.833, t[1]);
            double Dhuhr = this.computeMidDay(t[2]);
            double Asr = this.computeAsr(1 + this.asrJuristic, t[3]);
            double Sunset = this.computeTime(0.833, t[4]); ;
            double Maghrib = this.computeTime(this.methodParams[this.calcMethod][2], t[5]);
            double Isha = this.computeTime(this.methodParams[this.calcMethod][4], t[6]);

            return new double[] { Fajr, Sunrise, Dhuhr, Asr, Sunset, Maghrib, Isha };
        }


        // compute prayer times at given julian date
        private string[] computeDayTimes()
        {
            double[] times = new double[] { 5, 6, 12, 13, 18, 18, 18 }; //default times

            for (int i = 1; i <= this.numIterations; i++)
                times = this.computeTimes(times);

            times = this.adjustTimes(times);
            return this.adjustTimesFormat(times);
        }


        // adjust times in a prayer time array
        private double[] adjustTimes(double[] times)
        {
            for (int i = 0; i < 7; i++)
                times[i] += this.timeZone - this.lng / 15;
            times[2] += this.dhuhrMinutes / 60; //Dhuhr
            if (this.methodParams[this.calcMethod][1] == 1) // Maghrib
                times[5] = times[4] + this.methodParams[this.calcMethod][2] / 60;
            if (this.methodParams[this.calcMethod][3] == 1) // Isha
                times[6] = times[5] + this.methodParams[this.calcMethod][4] / 60;

            if (this.adjustHighLats != this.None)
                times = this.adjustHighLatTimes(times);
            return times;
        }


        // convert times array to given time format
        private string[] adjustTimesFormat(double[] times)
        {
            string[] returnData = new string[times.Length];

            if (this.timeFormat == this.Float)
                return returnData;
            for (int i = 0; i < 7; i++)
                if (this.timeFormat == this.Time12)
                    returnData[i] = floatToTime12(times[i]);
                else if (this.timeFormat == this.Time12NS)
                    returnData[i] = floatToTime12(times[i], true);
                else
                    returnData[i] = floatToTime24(times[i]);
            return returnData;
        }


        // adjust Fajr, Isha and Maghrib for locations in higher latitudes
        private double[] adjustHighLatTimes(double[] times)
        {
            double nightTime = this.timeDiff(times[4], times[1]); // sunset to sunrise

            // Adjust Fajr
            double FajrDiff = this.nightPortion(methodParams[this.calcMethod][0]) * nightTime;
            if (double.IsNaN(times[0]) || this.timeDiff(times[0], times[1]) > FajrDiff)
                times[0] = times[1] - FajrDiff;

            // Adjust Isha
            double IshaAngle = (this.methodParams[this.calcMethod][3] == 0) ? this.methodParams[this.calcMethod][4] : 18;
            double IshaDiff = this.nightPortion(IshaAngle) * nightTime;
            if (double.IsNaN(times[6]) || this.timeDiff(times[4], times[6]) > IshaDiff)
                times[6] = times[4] + IshaDiff;

            // Adjust Maghrib
            double MaghribAngle = (this.methodParams[this.calcMethod][1] == 0) ? this.methodParams[this.calcMethod][2] : 4;
            double MaghribDiff = this.nightPortion(MaghribAngle) * nightTime;
            if (double.IsNaN(times[5]) || this.timeDiff(times[4], times[5]) > MaghribDiff)
                times[5] = times[4] + MaghribDiff;

            return times;
        }


        // the night portion used for adjusting times in higher latitudes
        private double nightPortion(double angle)
        {
            if (this.adjustHighLats == this.AngleBased)
                return 1 / 60 * angle;
            if (this.adjustHighLats == this.MidNight)
                return 1 / 2d;
            if (this.adjustHighLats == this.OneSeventh)
                return 1 / 7d;

            return 0;
        }


        // convert hours to day portions 
        private double[] dayPortion(double[] times)
        {
            for (int i = 0; i < 7; i++)
                times[i] /= 24;
            return times;
        }


        #endregion

        #region Misc Functions
        // compute the difference between two times 
        private double timeDiff(int time1, int time2)
        {
            return this.fixhour(time2 - time1);
        }
        private double timeDiff(double time1, double time2)
        {
            return this.fixhour(time2 - time1);
        }

        // add a leading 0 if necessary
        private string twoDigitsFormat(int num)
        {
            return (num < 10) ? "0" + num.ToString() : num.ToString();
        }
        #endregion

        #region Julian Date Functions
        // calculate julian date from a calendar date
        private double julianDate(int year, int month, int day)
        {
            double A = Math.Floor((double)(year / 100));
            double B = Math.Floor(A / 4);
            double C = 2 - A + B;
            double Dd = day;
            double Ee = Math.Floor(365.25 * (year + 4716));
            double F = Math.Floor(30.6001 * (month + 1));

            double JD = C + Dd + Ee + F - 1524.5;

            return JD;

            //if (month <= 2)
            //{
            //    year -= 1;
            //    month += 12;
            //}
            //double A = Math.Floor(year / 100D);
            //double B = 2 - A + Math.Floor(A / 4);

            //double JD = Math.Floor(365.25 * (year + 4716)) + Math.Floor(30.6001 * (month + 1)) + day + B - 1524.5;
            //return JD;
        }


        // convert a calendar date to julian date (second method)
        private double calcJD(int year, int month, int day)
        {
            double J1970 = 2440588.0;
            DateTime date = new DateTime(year, month - 1, day);
            TimeSpan TS = new TimeSpan(year, month - 1, day);
            double ms = TS.TotalMilliseconds;   // # of milliseconds since midnight Jan 1, 1970
            double days = Math.Floor((double)ms / (1000 * 60 * 60 * 24));
            return J1970 + days - 0.5;
        }
        #endregion

        #region Time-Zone Functions
        // compute local time-zone for a specific date
        private int getTimeZone(DateTime date)
        {
            DateTime localDate = new DateTime(date.Year, date.Month, date.Day, 0, 0, 0, 0);
            string GMTString = localDate.ToString("R");
            DateTime GMTDate = DateTime.Parse(GMTString.Substring(0, GMTString.LastIndexOf(" ") - 1));
            int hoursDiff = (localDate - GMTDate).Hours / (1000 * 60 * 60);
            return hoursDiff;
        }


        // compute base time-zone of the system
        private int getBaseTimeZone()
        {
            return this.getTimeZone(new DateTime(2000, 0, 15));
        }


        // detect daylight saving in a given date
        private bool detectDaylightSaving(DateTime date)
        {
            return this.getTimeZone(date) != this.getBaseTimeZone();
        }


        // return effective timezone for a given date
        private int effectiveTimeZone(int year, int month, int day, int? timeZone)
        {
            if (timeZone == null)
                timeZone = this.getTimeZone(new DateTime(year, month - 1, day));
            return 1 * timeZone.Value;
        }

        #endregion

        #region Trigonometric Functions
        // degree sin
        private double dsin(double d)
        {
            return Math.Sin(this.dtr(d));
        }

        // degree cos
        private double dcos(double d)
        {
            return Math.Cos(this.dtr(d));
        }

        // degree tan
        private double dtan(double d)
        {
            return Math.Tan(this.dtr(d));
        }

        // degree arcsin
        private double darcsin(double x)
        {
            return this.rtd(Math.Asin(x));
        }

        // degree arccos
        private double darccos(double x)
        {
            return this.rtd(Math.Acos(x));
        }

        // degree arctan
        private double darctan(double x)
        {
            return this.rtd(Math.Atan(x));
        }

        // degree arctan2
        private double darctan2(double y, double x)
        {
            return this.rtd(Math.Atan2(y, x));
        }

        // degree arccot
        private double darccot(double x)
        {
            return this.rtd(Math.Atan(1 / x));
        }

        // degree to radian
        private double dtr(double d)
        {
            return (d * Math.PI) / 180.0;
        }

        // radian to degree
        private double rtd(double r)
        {
            return (r * 180.0) / Math.PI;
        }

        // range reduce angle in degrees.
        private double fixangle(double a)
        {
            a = a - 360.0 * (Math.Floor(a / 360.0));
            a = a < 0 ? a + 360.0 : a;
            return a;
        }

        // range reduce hours to 0..23
        private double fixhour(double a)
        {
            a = a - 24.0 * (Math.Floor(a / 24.0));
            a = a < 0 ? a + 24.0 : a;
            return a;
        }

        #endregion
    }


}
