#include <SD.h>
#include <Wire.h>
#include "RTClib.h"

// A simple data logger for the Arduino an analog pin

// how many milliseconds between grabbing data and logging it. 1000 ms is once a second
#define LOG_INTERVAL  100 // mills between entries (reduce to take more/faster data)

// how many milliseconds before writing the logged data permanently to disk
// set it to the LOG_INTERVAL to write each time (safest)
// set it to 10*LOG_INTERVAL to write all data every 10 datareads, you could lose up to 
// the last 10 reads if power is lost but it uses less power and is much faster!
#define SYNC_INTERVAL 30*LOG_INTERVAL // mills between calls to flush() - to write data to the card
unsigned long syncTime = 0; // time of last sync()

// Initialize the serial port and echo data and error messsages
#define ECHO_TO_SERIAL   0

// the digital pins that connect to the LEDs
#define redLEDpin 8
#define greenLEDpin 7
#define RangeTriggerPin 9
// The analog pins that connect to the sensors
#define RangePin 0           // analog 0
const byte Rangeis[7] = {0,1,2,3,6,5,4};

RTC_DS1307 RTC; // define the Real Time Clock object

// for the data logging shield, we use digital pin 10 for the SD cs line
const byte chipSelect = 10;

// the logging file
File logfile;

void error(char *str)
{
  // red LED indicates error
  digitalWrite(redLEDpin, HIGH);
  // Describe the error
#if ECHO_TO_SERIAL
  Serial.print("error: ");
  Serial.println(str);

#endif  //ECHO_TO_SERIAL

  logfile.print("error: ");
  logfile.println(str);  

  while(1);
}

void showrange(byte meters)
{
  for (byte i = 0; i < 7; i++) {
    if(i<meters){
     digitalWrite(Rangeis[i],HIGH);
    }
    else{
      digitalWrite(Rangeis[i],LOW);
    }
  }
}

void setup(void)
{
#if ECHO_TO_SERIAL
  // Initialize the serial port
  Serial.begin(9600);
#endif //ECHO_TO_SERIAL
  
  // use debugging LEDs
  pinMode(redLEDpin, OUTPUT);
  pinMode(greenLEDpin, OUTPUT);

  // Set ranging pin
  pinMode(RangeTriggerPin, OUTPUT);
  digitalWrite(RangeTriggerPin, HIGH);
  delay(100);

  for (byte i = 0; i < 7; i++) {
    pinMode(Rangeis[i],OUTPUT);
  }
  
  // connect to RTC
  Wire.begin();  
  // initialize the SD card make sure that the default chip select pin is set to output, even if you don't use it:
  pinMode(10, OUTPUT);  
  // see if the card is present and can be initialized:
  if (!SD.begin(chipSelect)) {
     error("Card failed, or not present");
  }
#if ECHO_TO_SERIAL
  Serial.println("card initialized.");
#endif  //ECHO_TO_SERIAL  
  // create a new file
  char filename[] = "LOGGER00.CSV";
  for (uint8_t i = 0; i < 100; i++) {
    filename[6] = i/10 + '0';
    filename[7] = i%10 + '0';
    if (! SD.exists(filename)) {
      // only open a new file if it doesn't exist
      logfile = SD.open(filename, FILE_WRITE); 
      break;  // leave the loop!
    }
  }
  if (! logfile) {
    error("couldnt create file");
  }
  // Connect to the RTC
  if (!RTC.begin()) {
    error("RTC failed");
  }
#if ECHO_TO_SERIAL
  Serial.println("RTC initialized.");
#endif  //ECHO_TO_SERIAL
#if ECHO_TO_SERIAL  
  Serial.print("Logging to: ");
  Serial.println(filename);
#endif // ECHO_TO_SERIAL

  logfile.println("millis,stamp,datetime,Range(m),RawRange");    
#if ECHO_TO_SERIAL
  Serial.println("millis,stamp,datetime,Range(m),RawRange");    
#endif //ECHO_TO_SERIAL

}

void loop(void)
{
  DateTime now;

  // delay for the amount of time we want between readings
  delay((LOG_INTERVAL -1) - (millis() % LOG_INTERVAL));
  
  // Notify the viewer that the device is reding
  digitalWrite(greenLEDpin, HIGH);
  // milliseconds since starting
  unsigned long m = millis();

  // fetch the time
  now = RTC.now();

  analogRead(RangePin);
  delay(10);
  int RangeReading = analogRead(RangePin);
  digitalWrite(greenLEDpin, LOW);
    
  // converting that reading to SI units
  byte RangeDistance = RangeReading/100; // Range in meters
  showrange(RangeDistance);
  
  // log data
  logfile.print(m);           // milliseconds since start
  logfile.print(", ");    
  logfile.print(now.unixtime()); // seconds since 1/1/1970
  logfile.print(", ");
  logfile.print('"');
  logfile.print(now.year(), DEC);
  logfile.print("/");
  logfile.print(now.month(), DEC);
  logfile.print("/");
  logfile.print(now.day(), DEC);
  logfile.print(" ");
  logfile.print(now.hour(), DEC);
  logfile.print(":");
  logfile.print(now.minute(), DEC);
  logfile.print(":");
  logfile.print(now.second(), DEC);
  logfile.print('"');
  logfile.print(", ");    
  logfile.print(RangeDistance); //Record distance
  logfile.print(", ");    
  logfile.print(RangeReading);
  logfile.println();
#if ECHO_TO_SERIAL
  Serial.print(m);         // milliseconds since start
  Serial.print(", ");  
  Serial.print(now.unixtime()); // seconds since 1/1/1970
  Serial.print(", ");
  Serial.print('"');
  Serial.print(now.year(), DEC);
  Serial.print("/");
  Serial.print(now.month(), DEC);
  Serial.print("/");
  Serial.print(now.day(), DEC);
  Serial.print(" ");
  Serial.print(now.hour(), DEC);
  Serial.print(":");
  Serial.print(now.minute(), DEC);
  Serial.print(":");
  Serial.print(now.second(), DEC);
  Serial.print('"');
  Serial.print(", ");    
  Serial.print(RangeDistance); //Record Distance
  Serial.print(", ");    
  Serial.print(RangeReading); 
  Serial.println();
#endif // ECHO_TO_SERIAL

  // Now we write data to disk! Don't sync too often - requires 2048 bytes of I/O to SD card
  // which uses a bunch of power and takes time
  if ((millis() - syncTime) < SYNC_INTERVAL) return;
  syncTime = millis();
  
  // blink both LEDs to show we are syncing data to the card & updating FAT!
  digitalWrite(greenLEDpin, HIGH);
  digitalWrite(redLEDpin, HIGH);
  logfile.flush();
  digitalWrite(redLEDpin, LOW);
  digitalWrite(greenLEDpin, LOW);
  
}

