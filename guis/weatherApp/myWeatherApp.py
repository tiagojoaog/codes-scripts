#!/usr/bin/env python3

import tkinter as tk
import os, sys
from PIL import ImageTk, Image
import requests
import json
import datetime
from functools import partial


class Weather():
    def __init__(self):
        self.win       = {}
        self.widgets   = {}
        self.data      = {}
        self.countries = ["AF", "AX", "AL", "DZ", "AS", "AD", "AO", "AI", "AQ", "AG", "AR", 
        "AM", "AW", "AU", "AT", "AZ", "BS", "BH", "BD", "BB", "BY", "BE", "BZ", "BJ", "BM", 
        "BT", "BO", "BQ", "BA", "BW", "BV", "BR", "IO", "BN", "BG", "BF", "BI", "KH", "CM", 
        "CA", "CV", "KY", "CF", "TD", "CL", "CN", "CX", "CC", "CO", "KM", "CG", "CD", "CK", 
        "CR", "CI", "HR", "CU", "CW", "CY", "CZ", "DK", "DJ", "DM", "DO", "EC", "EG", "SV", 
        "GQ", "ER", "EE", "ET", "FK", "FO", "FJ", "FI", "FR", "GF", "PF", "TF", "GA", "GM", 
        "GE", "DE", "GH", "GI", "GR", "GL", "GD", "GP", "GU", "GT", "GG", "GN", "GW", "GY", 
        "HT", "HM", "VA", "HN", "HK", "HU", "IS", "IN", "ID", "IR", "IQ", "IE", "IM", "IL", 
        "IT", "JM", "JP", "JE", "JO", "KZ", "KE", "KI", "KP", "KR", "KW", "KG", "LA", "LV", 
        "LB", "LS", "LR", "LY", "LI", "LT", "LU", "MO", "MK", "MG", "MW", "MY", "MV", "ML", 
        "MT", "MH", "MQ", "MR", "MU", "YT", "MX", "FM", "MD", "MC", "MN", "ME", "MS", "MA", 
        "MZ", "MM", "NA", "NR", "NP", "NL", "NC", "NZ", "NI", "NE", "NG", "NU", "NF", "MP", 
        "NO", "OM", "PK", "PW", "PS", "PA", "PG", "PY", "PE", "PH", "PN", "PL", "PT", "PR", 
        "QA", "RE", "RO", "RU", "RW", "BL", "SH", "KN", "LC", "MF", "PM", "VC", "WS", "SM",
        "ST", "SA", "SN", "RS", "SC", "SL", "SG", "SX", "SK", "SI", "SB", "SO", "ZA", "GS", 
        "SS", "ES", "LK", "SD", "SR", "SJ", "SZ", "SE", "CH", "SY", "TW", "TJ", "TZ", "TH", 
        "TL", "TG", "TK", "TO", "TT", "TN", "TR", "TM", "TC", "TV", "UG", "UA", "AE", "GB", 
        "US", "UM", "UY", "UZ", "VU", "VE", "VN", "VG", "VI", "WF", "EH", "YE", "ZM", "ZW"] 
    

        self.languages = {
                        "en" : "English",
                        "ar" : "Arabic",
                        "az" : "Azerbaijani",
                        "be" : "Belarusian",
                        "bg" : "Bulgarian",
                        "bs" : "Bosnian",
                        "ca" : "Catalan",
                        "cz" : "Czech",
                        "da" : "Danish",
                        "de" : "German",
                        "fi" : "Finnish",
                        "fr" : "French",
                        "el" : "Greek",
                        "et" : "Estonian",
                        "hr" : "Croation",
                        "hu" : "Hungarian",
                        "id" : "Indonesian",
                        "it" : "Italian",
                        "is" : "Icelandic",
                        "kw" : "Cornish",
                        "nb" : "Norwegian Bokmål",
                        "nl" : "Dutch",
                        "pl" : "Polish",
                        "pt" : "Portuguese",
                        "ro" : "Romanian",
                        "ru" : "Russian",
                        "sk" : "Slovak",
                        "sl" : "Slovenian",
                        "sr" : "Serbian",
                        "sv" : "Swedish",
                        "tr" : "Turkish",
                        "uk" : "Ukrainian",
                        "zh" : "Chinese (Simplified)",
                        "zh-tw" : "Chinese (Traditional)"
                }



        self.win["root"] = tk.Tk()
        self.win["root"].title("WeatherApp")
        self.win["root"].geometry("450x100")

        self.widgets["postal_code"] = tk.Entry(self.win["root"])
        self.widgets["postal_code"].insert(0, "Write postal code here")

        self.widgets["countries"] = tk.Entry(self.win['root'])
        self.widgets["countries"].insert(0, "Write country code here.")

        self.widgets["latitude"] = tk.Entry(self.win['root'])
        self.widgets["latitude"].insert(0, "Write latitude here.")

        self.widgets["longitude"] = tk.Entry(self.win['root'])
        self.widgets["longitude"].insert(0, "Write longitude here.")

        self.widgets["retrieve_button"] = tk.Button(self.win["root"], text="Retrieve", command=self.retrieve)
        
        self.widgets["language_variable"] = tk.StringVar(self.win['root'])
        self.widgets["language"] = tk.OptionMenu(self.win['root'], self.widgets["language_variable"], *self.languages.values())
        self.widgets["language_variable"].set(self.languages["en"])


        self.widgets["postal_code"].grid(row=0, column=0)
        self.widgets["language"].grid(row=0, column=2)
        self.widgets["countries"].grid(row=0, column=1)
        self.widgets["latitude"].grid(row=1, column=0)
        self.widgets["longitude"].grid(row=1, column=1)
        self.widgets["retrieve_button"].grid(row=2, column=0,columnspan=4)


    def search_key(self,_dict,name):
        for index,value in _dict.items():
            if name == value:
                return index


    def color_range(self,variable1,variable2):
        if min(variable1,variable2) <0:
            return "Blue4"
        elif max(variable1, variable2) >35:
            return "Maroon"
        else:
            return "Green"

    def read_data(self,day,event=None):
        self.clear()

        # Current_data
        ultraviolet=int(self.api_current["data"][day]["uv"])
        if ultraviolet <= 2:
            ultraviolet=f"{ultraviolet} (low)"
            ucolor="Green"
        elif ultraviolet <= 5:
            ultraviolet=f"{ultraviolet} (moderate)"
            ucolor="Darkorange1"
        elif ultraviolet <= 7:
            ultraviolet=f"{ultraviolet} (High)"
            ucolor="Red"
        elif ultraviolet <= 10:
            ultraviolet=f"{ultraviolet} (Very high)"
            ucolor="Maroon"
        else:
            ultraviolet=f"{ultraviolet} (Extreme)"
            ucolor="Purple"


        precipitation = round(float(self.api_current["data"][day]["pop"]),1) # %
        temperature   = round(float(self.api_current["data"][day]["temp"]),1)   # °C
        max_temp      = round(float(self.api_current["data"][day]["max_temp"]),1)
        min_temp      = round(float(self.api_current["data"][day]["min_temp"]),1)
        feels_like_max= round(float(self.api_current["data"][day]["app_max_temp"]),1)   # °C
        feels_like_min= round(float(self.api_current["data"][day]["app_min_temp"]),1)
        wind_maxspeed = round(float(self.api_current["data"][day-1]["wind_gust_spd"]),2)   #km/h
        wind_dir      = self.api_current["data"][day]["wind_cdir_full"]
        clouds        = round(float(self.api_current["data"][day]["clouds"]),1) # %
        visibility    = round(float(self.api_current["data"][day]["vis"]),1) # Km
        humidity      = self.api_current["data"][day]["rh"] # %
        snow          = self.api_current["data"][day]["snow_depth"] # mm
        pressure      = self.api_current["data"][day]["pres"]  #mbar
        date          = self.api_current["data"][day]["valid_date"]
        details       = self.api_current["data"][day]["weather"]["description"]

        tcolor  =self.color_range(min_temp,max_temp)
        tfcolor =self.color_range(feels_like_min,feels_like_max)
        
        self.widgets[f"label7"]=tk.Label(self.win["weather"], text=f"Date: {date} (Local) \n", justify=tk.LEFT, font="Helvetica 12")
        self.widgets[f"label8"]=tk.Label(self.win["weather"], text=f"Sun and Temperature:", justify=tk.LEFT, font="Helvetica 12")
        self.widgets[f"label9"]=tk.Label(self.win["weather"], text=f"  Temperature  {min_temp} (min.) {max_temp} (max.) {temperature} (avg.) °C", justify=tk.LEFT, fg=tcolor)
        self.widgets[f"label10"]=tk.Label(self.win["weather"], text=f"  Feels like  {feels_like_min} (min.) {feels_like_max} (max.) °C", justify=tk.LEFT, fg=tfcolor)
        self.widgets[f"label11"]=tk.Label(self.win["weather"], text=f"  Ultraviolet index  {ultraviolet}\n", justify=tk.LEFT, fg=ucolor)
        self.widgets[f"label12"]=tk.Label(self.win["weather"], text=f"Weather Phenomena:", justify=tk.LEFT, font="Helvetica 12")
        self.widgets[f"label13"]=tk.Label(self.win["weather"], text=f"  Precipitation chance  {precipitation} %", justify=tk.LEFT)
        self.widgets[f"label14"]=tk.Label(self.win["weather"], text=f"  Snow depth  {snow} mm\n", justify=tk.LEFT)
        self.widgets[f"label15"]=tk.Label(self.win["weather"], text=f"Clouds and Air:", justify=tk.LEFT, font="Helvetica 12")
        self.widgets[f"label16"]=tk.Label(self.win["weather"], text=f"  Humidity  {humidity} %", justify=tk.LEFT)
        self.widgets[f"label17"]=tk.Label(self.win["weather"], text=f"  Wind  {wind_maxspeed} Km/h ({wind_dir})", justify=tk.LEFT)
        self.widgets[f"label18"]=tk.Label(self.win["weather"], text=f"  Cloud coverage  {clouds} %", justify=tk.LEFT)
        self.widgets[f"label19"]=tk.Label(self.win["weather"], text=f"  Visibility  {visibility} Km", justify=tk.LEFT)
        self.widgets[f"label20"]=tk.Label(self.win["weather"], text=f"  Pressure  {pressure} mbar\n", justify=tk.LEFT)
        self.widgets[f"label21"]=tk.Label(self.win["weather"], text=f"Details: {details}", justify=tk.LEFT, font="Helvetica 12")

        for i in range(15):
            i+=7
            self.widgets[f"label{i}"].grid(row=i,column=0,columnspan=20, sticky=tk.W)


    def clear(self):
        for i in range(30):
            i+=7
            label = f"label{i}"
            if label in self.widgets:
                try:
                    self.widgets[label].config(text="")
                except:
                    pass
                
    def big_spaces(self):
        for i in range(10):
            self.widgets[f"label{i}"]=tk.Label(self.win["weather"], text="")
            self.widgets[f"label{i}"].grid(row=6+i,column=0)
   
    def weekday(self,day,i):
        while day>6:
            day-=7

        if i==0:
            return "Today"
        elif day==0:
            return "Mon"
        elif day==1:
            return "Tue"
        elif day==2:
            return "Wed"
        elif day==3:
            return "Thu"
        elif day==4:
            return "Fri"
        elif day==5:
            return "Sat"
        elif day==6:
            return "Sun"

    def click(self,i=0, event=None):
        if "warnings" in self.widgets:
            if self.widgets["warnings"]:
                self.widgets["warnings"].destroy()

        alert=self.api_alert["alerts"]
        self.alert_message =f'''    ***ALERT***
    Title:  {self.api_alert["alerts"][i]["title"]}

    Description: {self.api_alert["alerts"][i]["description"]}

    Severity: {self.api_alert["alerts"][i]["severity"]}

    Locations: {str(self.api_alert["alerts"][i]["regions"])[1:-1]}    '''
        if i< len(alert)-1 and i> 0:
            statusr=tk.NORMAL
            statusl=tk.NORMAL
        elif i==0:
            statusl=tk.DISABLED
            statusr=tk.NORMAL
        else:
            statusr=tk.DISABLED
            statusl=tk.NORMAL

        
        alert_color="Red"
        self.widgets['warnings']=tk.Label(self.win["weather"], text=self.alert_message, fg=alert_color, bg="White",justify="center", borderwidth=3, relief="sunken")
        self.widgets['warnings'].grid(row=38,column=0,columnspan=20)
        
        tk.Label(self.win["weather"], text="").grid(row=39,column=0,columnspan=20)

        tk.Button(self.win["weather"], text="<", state=statusl, command=partial(self.click,i=i-1)).grid(row=40,column=3)
        tk.Button(self.win["weather"], text=">", state=statusr, command=partial(self.click,i=i+1)).grid(row=40,column=4)

    def retrieve(self,event=None):
        if "weather" in self.win:
            if self.win["weather"]:
                self.win["weather"].destroy()
                self.win["weather"]= None

    
        self.api_alert ={"alerts" :[]}

        self.win["weather"] = tk.Toplevel()
        self.win["weather"].title("WeatherApp")
        self.win["weather"].geometry("800x1000")

        try:
            float(self.widgets['latitude'].get())
            float(self.widgets['longitude'].get())
            coordinates=True
        except ValueError:
            coordinates=False


        country= str(self.widgets['countries'].get()).upper()
        try:
            lang = self.search_key(self.languages,self.widgets["language_variable"].get())
            
            if coordinates:
                api_request_current = requests.get(f"https://api.weatherbit.io/v2.0/forecast/daily?lat={self.widgets['latitude'].get()}&lon={self.widgets['longitude'].get()}&lang={lang}&key=131539a0d631498c8dde6333afa424b7")
            else:
                api_request_current = requests.get(f"https://api.weatherbit.io/v2.0/forecast/daily?postal_code={self.widgets['postal_code'].get()}&country={country}&lang={lang}&key=131539a0d631498c8dde6333afa424b7")

            
            text_color="Gray1"

            self.api_current      = json.loads(api_request_current.content)
                
            City          = self.api_current["city_name"]
            CountryCode   = self.api_current["country_code"]
            timezone      = self.api_current["timezone"]
            longitude     = self.api_current["lon"] #degrees
            latitude      = self.api_current["lat"] #degrees

            api_request_now       = requests.get(f"https://api.weatherbit.io/v2.0/current?lat={latitude}&lon={longitude}&lang={lang}&key=131539a0d631498c8dde6333afa424b7")
            self.api_now          = json.loads(api_request_now.content)

            current_temp  = round(float(self.api_now["data"][0]["temp"]),1)
            current_precip= self.api_now["data"][0]["precip"]
            current_descr = self.api_now["data"][0]["weather"]["description"]
            
            try:
                api_request_alert = requests.get(f"https://api.weatherbit.io/v2.0/alerts?lat={latitude}&lon={longitude}&key=131539a0d631498c8dde6333afa424b7")
                self.api_alert    = json.loads(api_request_alert.content)
                self.alert_message =f'''    ***ALERT*** 
    Title:  {self.api_alert["alerts"][0]["title"]}
    
    Description: {self.api_alert["alerts"][0]["description"]}
    
    Severity: {self.api_alert["alerts"][0]["severity"]}
    
    Locations: {str(self.api_alert["alerts"][0]["regions"])[1:-1]}    '''
                alert_color = "Red"
                self.win["weather"].geometry("800x900")
            except:
                self.alert_message = '''    ***ALERT*** 
    Nothing to declare.   '''
                alert_color   = "Gray1"

            self.widgets["countries"].delete(0, tk.END)
            self.widgets["latitude"].delete(0, tk.END)
            self.widgets["longitude"].delete(0, tk.END)

            self.widgets["countries"].insert(0, CountryCode)
            self.widgets["latitude"].insert(0, latitude)
            self.widgets["longitude"].insert(0, longitude)
             
            self.header1 = tk.Label(self.win["weather"], text=f"Location: {City} ({CountryCode} lat {round(float(latitude),4)} long {round(float(longitude),4)})", font="Helvetica 18 bold").grid(row=0,column=0,columnspan=20)
            self.header2 = tk.Label(self.win["weather"], text=f"Timezone: {timezone}", font="Helvetica 16").grid(row=1,column=0,columnspan=20)
            self.header3 = tk.Label(self.win["weather"], text=f"{current_temp} °C  Chance of rain {current_precip} % ({current_descr})", font="Helvetica 16").grid(row=2,column=0,columnspan=20)

            tk.Label(self.win["weather"], text="").grid(row=3,column=0)

            self.make_buttons = [tk.Button(self.win["weather"], text=self.weekday(datetime.datetime.today().weekday()+i,i), width=7, command=partial(self.read_data,i)) for i in range(16)]
            
            for index,item in enumerate(self.make_buttons):
                if index>=8:
                    row=5
                    column=index-8

                else:
                    row=4
                    column=index

                item.grid(row=row,column=column)
            
            self.big_spaces()

            self.message      = f"Weather information retrieved from weatherbit.io"
            
            self.widgets["latitude"].delete(0, tk.END)
            self.widgets["longitude"].delete(0, tk.END)

            self.widgets["latitude"].insert(0, "Write latitude here.")
            self.widgets["longitude"].insert(0, "Write longitude here.")
         

        except Exception as msgerror:
            text_color="Red"
            self.alert_message =''
            alert_color   = "Gray1"
            if country not in self.countries and not coordinates:
                self.message="ERROR: Country code does not exist. Make sure you write the correct 2-letter code."
            
            else:
                self.message = f'''
                ERROR: Could not retrieve data. Postal Code, Country or lat/long may not exist in the database.
                    If Postal code is stated, country has to be as well, otherwise, latitude and longitude will suffice.
                    More information about this error --> ({msgerror})
                '''
        self.widgets["warnings"]=tk.Label(self.win["weather"], text=self.alert_message, fg=alert_color, bg="White",justify="center", borderwidth=3, relief="sunken")
        self.widgets["warnings"].grid(row=38,column=0,columnspan=20)
        tk.Label(self.win["weather"], text="").grid(row=39,column=0,columnspan=20)
        
        if len(self.api_alert["alerts"])<2:
            state= tk.DISABLED
        else:
            state= tk.NORMAL

        tk.Button(self.win["weather"], text="<", state=tk.DISABLED).grid(row=40,column=3)
        tk.Button(self.win["weather"], text=">", state=state, command=partial(self.click,i=1)).grid(row=40,column=4)
        tk.Label(self.win["weather"], text="").grid(row=41,column=0,columnspan=20)
        tk.Label(self.win["weather"], text=self.message, fg=text_color).grid(row=42,column=0,columnspan=20)

def main():
    
    WeatherObj = Weather()

    WeatherObj.win["root"].mainloop()



if __name__ == "__main__":
    main()
