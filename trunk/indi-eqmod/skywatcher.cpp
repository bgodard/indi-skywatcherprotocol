/* Copyright 2012 Geehalel (geehalel AT gmail DOT com) */
/* This file is part of the Skywatcher Protocol INDI driver.

    The Skywatcher Protocol INDI driver is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Skywatcher Protocol INDI driver is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Skywatcher Protocol INDI driver.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <termios.h>
#include <memory>

#include <libnova.h>

/* Normally unused-only for debug */
#include <indidevapi.h>
#include <indicom.h>

#include "eqmoderror.h"
#include "skywatcher.h"

Skywatcher::Skywatcher(void)
{
  fd=-1;

}

Skywatcher::~Skywatcher(void)
{
    Disconnect();
}


/* API */

bool Skywatcher::Connect(char *port)  throw (EQModError)
{   
  int err_code = 0;
  unsigned long tmpMCVersion=0;
  if ((err_code=tty_connect(port, 9600, 8, 0, 1, &fd)) != TTY_OK)
    {
      char ttyerrormsg[ERROR_MSG_LENGTH];
      tty_error_msg(err_code, ttyerrormsg, ERROR_MSG_LENGTH);
      throw EQModError(EQModError::ErrDisconnect, "Error connecting to port %s: %s", 
		       port, ttyerrormsg);
      return false;
    }
    
  dispatch_command(InquireMotorBoardVersion, Axis1, NULL);
  read_eqmod();
  tmpMCVersion=Revu24str2long(response+1);
  MCVersion = ((tmpMCVersion & 0xFF) << 16) | ((tmpMCVersion & 0xFF00)) | ((tmpMCVersion & 0xFF0000) >> 16);
  MountCode=MCVersion & 0xFF;
  /* Check supported mounts here */
  if ((MountCode == 0x80) || (MountCode == 0x81) /*|| (MountCode == 0x82)*/ || (MountCode == 0x90)) {
    
    throw EQModError(EQModError::ErrDisconnect, "Mount not supported: mount code 0x%x (0x80=GT, 0x81=MF, 0x82=114GT, 0x90=DOB)", 
		     MountCode);
    return false;
  }

  return true;
}


bool Skywatcher::Disconnect() throw (EQModError)
{
  if (fd < 0) return true;
  StopMotor(Axis1);
  StopMotor(Axis2);
  // Deactivate motor (for geehalel mount only)
  if (MountCode == 0xF0) {
    dispatch_command(Deactivate, Axis1, NULL);
    read_eqmod();
  }
  tty_disconnect(fd);
  fd=-1;
  return true;
}

unsigned long Skywatcher::GetRAEncoder()  throw (EQModError)
{
  // Axis Position
  dispatch_command(GetAxisPosition, Axis1, NULL);
  read_eqmod();
  RAStep=Revu24str2long(response+1);
  gettimeofday (&lastreadmotorposition[Axis1], NULL);
  return RAStep;
}

unsigned long Skywatcher::GetDEEncoder()  throw (EQModError)
{  
  // Axis Position
  dispatch_command(GetAxisPosition, Axis2, NULL);
  read_eqmod();
 
  DEStep=Revu24str2long(response+1);
  gettimeofday (&lastreadmotorposition[Axis2], NULL);
  return DEStep;
}

unsigned long Skywatcher::GetRAEncoderZero()
{
  return RAStepInit;
}

unsigned long Skywatcher::GetRAEncoderTotal()
{
  return RASteps360;
}

unsigned long Skywatcher::GetDEEncoderZero()
{
  return DEStepInit;
}

unsigned long Skywatcher::GetDEEncoderTotal()
{
  return DESteps360;
}

unsigned long Skywatcher::GetRAPeriod() throw (EQModError)
{
  return RAPeriod;
}

unsigned long Skywatcher::GetDEPeriod() throw (EQModError)
{
  return DEPeriod;
}

void Skywatcher::GetRAMotorStatus(ILightVectorProperty *motorLP) throw (EQModError)
{

  ReadMotorStatus(Axis1);
  if (!RAInitialized) {
    IUFindLight(motorLP, "RAInitialized")->s=IPS_ALERT;
    IUFindLight(motorLP, "RARunning")->s=IPS_IDLE;
    IUFindLight(motorLP, "RAGoto")->s=IPS_IDLE;
    IUFindLight(motorLP, "RAForward")->s=IPS_IDLE;
    IUFindLight(motorLP, "RAHighspeed")->s=IPS_IDLE;
  } else {
    IUFindLight(motorLP, "RAInitialized")->s=IPS_OK;
    IUFindLight(motorLP, "RARunning")->s=(RARunning?IPS_OK:IPS_BUSY);
    IUFindLight(motorLP, "RAGoto")->s=((RAStatus.slewmode==GOTO)?IPS_OK:IPS_BUSY);
    IUFindLight(motorLP, "RAForward")->s=((RAStatus.direction==FORWARD)?IPS_OK:IPS_BUSY);
    IUFindLight(motorLP, "RAHighspeed")->s=((RAStatus.speedmode==HIGHSPEED)?IPS_OK:IPS_BUSY);
  }
}

void Skywatcher::GetDEMotorStatus(ILightVectorProperty *motorLP) throw (EQModError)
{

  ReadMotorStatus(Axis2);
  if (!DEInitialized) {
    IUFindLight(motorLP, "DEInitialized")->s=IPS_ALERT;
    IUFindLight(motorLP, "DERunning")->s=IPS_IDLE;
    IUFindLight(motorLP, "DEGoto")->s=IPS_IDLE;
    IUFindLight(motorLP, "DEForward")->s=IPS_IDLE;
    IUFindLight(motorLP, "DEHighspeed")->s=IPS_IDLE;
  } else {
    IUFindLight(motorLP, "DEInitialized")->s=IPS_OK;
    IUFindLight(motorLP, "DERunning")->s=(DERunning?IPS_OK:IPS_BUSY);
    IUFindLight(motorLP, "DEGoto")->s=((DEStatus.slewmode==GOTO)?IPS_OK:IPS_BUSY);
    IUFindLight(motorLP, "DEForward")->s=((DEStatus.direction==FORWARD)?IPS_OK:IPS_BUSY);
    IUFindLight(motorLP, "DEHighspeed")->s=((DEStatus.speedmode==HIGHSPEED)?IPS_OK:IPS_BUSY);
  }
}


void Skywatcher::Init(ISwitchVectorProperty *parkSP) throw (EQModError) 
{
  ReadMotorStatus(Axis1);
  ReadMotorStatus(Axis2);
  if (!RAInitialized && !DEInitialized) {
    //Read initial stepper values 
    dispatch_command(GetAxisPosition, Axis1, NULL);
    read_eqmod();
    RAStepInit=Revu24str2long(response+1);
    dispatch_command(GetAxisPosition, Axis2, NULL);
    read_eqmod();
    DEStepInit=Revu24str2long(response+1);
    if (parkSP->sp[0].s==ISS_ON) { 
      //TODO get Park position, set corresponding encoder values, mark mount as parked 
      IDMessage("Skywatcher driver", "Mount in park position");
    } else {
      //mount is supposed to be in the home position (pointing Celestial Pole)
      IDMessage("Skywatcher driver", "Mount in home position");
      RAStepHome=RAStepInit;
      DEStepHome=DEStepInit + (DESteps360 / 4);
      char cmdarg[7];
      cmdarg[6]='\0';
      long2Revu24str(DEStepHome, cmdarg);
      dispatch_command(SetAxisPosition, Axis2, cmdarg);
      read_eqmod();
      //TODO mark mount as unparked
    }
    // Energize motors
    dispatch_command(Initialize, Axis1, NULL);
    read_eqmod();
    dispatch_command(Initialize, Axis2, NULL);
    read_eqmod();

  } else { // Mount already initialized by another driver / driver instance
    // use default configuration && leave unchanged encoder values
    RAStepInit=0x800000;
    DEStepInit=0x800000;
    RAStepHome=RAStepInit;
    DEStepHome=DEStepInit + (DESteps360 / 4);
    IDMessage("Skywatcher driver", "Motors already initialized");
    // TODO mark mount as unparked
  }

}

void Skywatcher::InquireBoardVersion(ITextVectorProperty *boardTP) throw (EQModError)
{
  unsigned long tmpMCVersion=0;
  unsigned nprop=0;
  char *boardinfo[2];
  const char *boardinfopropnames[]={"MOUNT_TYPE", "MOTOR_CONTROLLER"};

  /*
  dispatch_command(InquireMotorBoardVersion, Axis1, NULL);
  read_eqmod();
  tmpMCVersion=Revu24str2long(response+1);
  MCVersion = ((tmpMCVersion & 0xFF) << 16) | ((tmpMCVersion & 0xFF00)) | ((tmpMCVersion & 0xFF0000) >> 16);
  MountCode=MCVersion & 0xFF;
  */
  nprop=2;
  //  strcpy(boardinfopropnames[0],"MOUNT_TYPE");
  boardinfo[0]=(char *) malloc(20*sizeof(char));
  switch(MountCode) {
  case 0x00: strcpy(boardinfo[0],"EQ6"); break;
  case 0x01: strcpy(boardinfo[0],"HEQ5"); break;
  case 0x02: strcpy(boardinfo[0],"EQ5"); break;
  case 0x03: strcpy(boardinfo[0],"EQ3"); break;
  case 0x80: strcpy(boardinfo[0],"GT"); break;
  case 0x81: strcpy(boardinfo[0],"MF"); break;
  case 0x82: strcpy(boardinfo[0],"114GT"); break;
  case 0x90: strcpy(boardinfo[0],"DOB"); break;
  case 0xF0: strcpy(boardinfo[0],"GEEHALEL"); break;
  default:  strcpy(boardinfo[0],"CUSTOM"); break;
  }
  
  boardinfo[1]=(char *) malloc(5);
  sprintf(boardinfo[1],"%04lx", (MCVersion >> 8));
  boardinfo[1][4]='\0';
  // should test this is ok
  IUUpdateText(boardTP, boardinfo, (char **)boardinfopropnames, nprop);
  IDSetText(boardTP,NULL);
  /* Check supported mounts here */
  /*if ((MountCode == 0x80) || (MountCode == 0x81) || (MountCode == 0x82) || (MountCode == 0x90)) {
    
    throw EQModError(EQModError::ErrDisconnect, "Mount not supported %s (mount code %d)", 
		     boardinfo[0], MountCode);
  }
  */
  free(boardinfo[0]); free(boardinfo[1]); 
}

void Skywatcher::InquireRAEncoderInfo(INumberVectorProperty *encoderNP) throw (EQModError) {
  double steppersvalues[3];
  const char *steppersnames[]={"RASteps360", "RAStepsWorm","RAHighspeedRatio"};
  // Steps per 360 degrees
  dispatch_command(InquireGridPerRevolution, Axis1, NULL);
  read_eqmod();
  RASteps360=Revu24str2long(response+1);
  steppersvalues[0]=(double)RASteps360;
  
  // Steps per Worm
  dispatch_command(InquireTimerInterruptFreq, Axis1, NULL);
  read_eqmod();
  RAStepsWorm=Revu24str2long(response+1);
 // There is a bug in the earlier version firmware(Before 2.00) of motor controller MC001.
  // Overwrite the GearRatio reported by the MC for 80GT mount and 114GT mount.
  if ((MCVersion & 0x0000FF) == 0x80)
    {
      RAStepsWorm = 0x162B97;           // for 80GT mount
    }
  if ((MCVersion & 0x0000FF) == 0x82)
    {
      RAStepsWorm = 0x205318;           // for 114GT mount
    }
  steppersvalues[1]=(double)RAStepsWorm;
  
  // Highspeed Ratio
  dispatch_command(InquireHighSpeedRatio, Axis1, NULL);
  read_eqmod();
  RAHighspeedRatio=Revu24str2long(response+1);
  
  steppersvalues[2]=(double)RAHighspeedRatio;
  // should test this is ok
  IUUpdateNumber(encoderNP, steppersvalues, (char **)steppersnames, 3);
  IDSetNumber(encoderNP, NULL);
}

void Skywatcher::InquireDEEncoderInfo(INumberVectorProperty *encoderNP) throw (EQModError) {
  double steppersvalues[3];
  const char *steppersnames[]={"DESteps360", "DEStepsWorm", "DEHighspeedRatio"};
  // Steps per 360 degrees
  dispatch_command(InquireGridPerRevolution, Axis2, NULL);
  read_eqmod();
  DESteps360=Revu24str2long(response+1);
  steppersvalues[0]=(double)DESteps360;
  
  // Steps per Worm
  dispatch_command(InquireTimerInterruptFreq, Axis2, NULL);
  read_eqmod();
  DEStepsWorm=Revu24str2long(response+1);
  // There is a bug in the earlier version firmware(Before 2.00) of motor controller MC001.
  // Overwrite the GearRatio reported by the MC for 80GT mount and 114GT mount.
  if ((MCVersion & 0x0000FF) == 0x80)
    {
      DEStepsWorm = 0x162B97;           // for 80GT mount
    }
  if ((MCVersion & 0x0000FF) == 0x82)
    {
      DEStepsWorm = 0x205318;           // for 114GT mount
    }

  steppersvalues[1]=(double)DEStepsWorm;

  // Highspeed Ratio
  dispatch_command(InquireHighSpeedRatio, Axis2, NULL);
  read_eqmod();
  DEHighspeedRatio=Revu24str2long(response+1);
  steppersvalues[2]=(double)DEHighspeedRatio;
  // should test this is ok
  IUUpdateNumber(encoderNP, steppersvalues, (char **)steppersnames, 3);
  IDSetNumber(encoderNP, NULL);
}

bool Skywatcher::IsRARunning() throw (EQModError)
{
  CheckMotorStatus(Axis1);
  return(RARunning);
}

bool Skywatcher::IsDERunning() throw (EQModError) 
{
  CheckMotorStatus(Axis2);
  return(DERunning);
}

void Skywatcher::ReadMotorStatus(SkywatcherAxis axis)  throw (EQModError)
{
  dispatch_command(GetAxisStatus, axis, NULL);
  read_eqmod();
  switch (axis) {
  case Axis1:
    RAInitialized=(response[3]&0x01);
    RARunning=(response[2]&0x01);
    if (response[1] & 0x01 ) RAStatus.slewmode=SLEW; else RAStatus.slewmode=GOTO;
    if (response[1] & 0x02 ) RAStatus.direction=BACKWARD; else RAStatus.direction=FORWARD;
    if (response[1] & 0x04 ) RAStatus.speedmode=HIGHSPEED; else RAStatus.speedmode=LOWSPEED;
    break;
  case Axis2:
    DEInitialized=(response[3]&0x01);
    DERunning=(response[2]&0x01);
    if (response[1] & 0x01 ) DEStatus.slewmode=SLEW; else DEStatus.slewmode=GOTO;
    if (response[1] & 0x02 ) DEStatus.direction=BACKWARD; else DEStatus.direction=FORWARD;
    if (response[1] & 0x04 ) DEStatus.speedmode=HIGHSPEED; else DEStatus.speedmode=LOWSPEED;
  default: ;
  }
  gettimeofday (&lastreadmotorstatus[axis], NULL);
}

void Skywatcher::SlewRA(double rate) throw (EQModError)
{
  double absrate=fabs(rate);
  unsigned long period=0;
  bool useHighspeed=false; 
  SkywatcherAxisStatus newstatus;

  
  if (RARunning && (RAStatus.slewmode==GOTO)) {
      throw  EQModError(EQModError::ErrInvalidCmd, "Can not slew while goto is in progress");
  } 

  if ((absrate < get_min_rate()) || (absrate > get_max_rate())) {
    throw  EQModError(EQModError::ErrInvalidParameter, "Speed rate out of limits: %.2fx Sidereal (min=%.2f, max=%.2f)", 
		      absrate, MIN_RATE, MAX_RATE);
  } 
  //if (MountCode != 0xF0) {
    if (absrate > SKYWATCHER_LOWSPEED_RATE) {
      absrate = absrate / RAHighspeedRatio;
      useHighspeed = true;
    }
    //}
  period=(long)(((SKYWATCHER_STELLAR_DAY * (double)RAStepsWorm) / (double)RASteps360) / absrate);
  if (rate >= 0.0) newstatus.direction = FORWARD; else newstatus.direction = BACKWARD;
  newstatus.slewmode=SLEW;
  if (useHighspeed) newstatus.speedmode = HIGHSPEED; else newstatus.speedmode = LOWSPEED;
  SetMotion(Axis1, newstatus);
  SetSpeed(Axis1, period);
  StartMotor(Axis1);
}

void Skywatcher::SlewDE(double rate) throw (EQModError)
{
  double absrate=fabs(rate);
  unsigned long period=0;
  bool useHighspeed=false; 
  SkywatcherAxisStatus newstatus;

  
  if (DERunning && (DEStatus.slewmode==GOTO)) {
      throw  EQModError(EQModError::ErrInvalidCmd, "Can not slew while goto is in progress");
  } 

  if ((absrate < get_min_rate()) || (absrate > get_max_rate())) {
    throw  EQModError(EQModError::ErrInvalidParameter, "Speed rate out of limits: %.2fx Sidereal (min=%.2f, max=%.2f)", 
		      absrate, MIN_RATE, MAX_RATE);
  } 
  //if (MountCode != 0xF0) {
    if (absrate > SKYWATCHER_LOWSPEED_RATE) {
      absrate = absrate / DEHighspeedRatio;
      useHighspeed = true;
    }
    //}
  period=(long)(((SKYWATCHER_STELLAR_DAY * (double)DEStepsWorm) / (double)DESteps360) / absrate);
  //IDLog("Slewing DE at %.2f %.2f %x %f\n", rate, absrate, period, (((SKYWATCHER_STELLAR_DAY * (double)RAStepsWorm) / (double)RASteps360) / absrate));
  if (rate >= 0.0) newstatus.direction = FORWARD; else newstatus.direction = BACKWARD;
  newstatus.slewmode=SLEW;
  if (useHighspeed) newstatus.speedmode = HIGHSPEED; else newstatus.speedmode = LOWSPEED;
  SetMotion(Axis2, newstatus);
  SetSpeed(Axis2, period);
  StartMotor(Axis2);
}

void Skywatcher::SlewTo(long deltaraencoder, long deltadeencoder)
{
  SkywatcherAxisStatus newstatus;
  bool useHighSpeed=false;
  unsigned long highperiod=10, lowperiod=18, lowspeedmargin = 20000, breaks=400;
  /* highperiod = RA 450X DE (+5) 200x, low period 32x */ 
  newstatus.slewmode=GOTO;
  if (deltaraencoder >= 0) newstatus.direction = FORWARD; else newstatus.direction = BACKWARD;
  if (deltaraencoder < 0) deltaraencoder = -deltaraencoder;
  if (deltaraencoder > lowspeedmargin) useHighSpeed = true; else useHighSpeed = false;
  if (useHighSpeed) newstatus.speedmode = HIGHSPEED; else newstatus.speedmode = LOWSPEED;
  if (deltaraencoder > 0) {
    SetMotion(Axis1, newstatus);
    if (useHighSpeed) SetSpeed(Axis1, highperiod); else SetSpeed(Axis1, lowperiod); 
    SetTarget(Axis1, deltaraencoder);
    if (useHighSpeed) breaks=((deltaraencoder > 3200) ? 3200 : deltaraencoder);
    else  breaks=((deltaraencoder > 200) ? 200 : deltaraencoder);
    SetTargetBreaks(Axis1, breaks);
    StartMotor(Axis1);
  }

  if (deltadeencoder >= 0) newstatus.direction = FORWARD; else newstatus.direction = BACKWARD;
  if (deltadeencoder < 0) deltadeencoder = -deltadeencoder;
  if (deltadeencoder > lowspeedmargin) useHighSpeed = true; else useHighSpeed = false;
  if (useHighSpeed) newstatus.speedmode = HIGHSPEED; else newstatus.speedmode = LOWSPEED;
  if (deltadeencoder > 0) {
    SetMotion(Axis2, newstatus);
    if (useHighSpeed) SetSpeed(Axis2, highperiod+5 ); else SetSpeed(Axis2, lowperiod); 
    SetTarget(Axis2, deltadeencoder);
    if (useHighSpeed) breaks=((deltadeencoder > 3200) ? 3200 : deltadeencoder);
    else  breaks=((deltadeencoder > 200) ? 200 : deltadeencoder);
    SetTargetBreaks(Axis2, breaks);
    StartMotor(Axis2);
  }

}

void  Skywatcher::SetRARate(double rate)  throw (EQModError)
{
  double absrate=fabs(rate);
  unsigned long period=0;
  bool useHighspeed=false; 
  SkywatcherAxisStatus newstatus;

  if ((absrate < get_min_rate()) || (absrate > get_max_rate())) {
    throw  EQModError(EQModError::ErrInvalidParameter, "Speed rate out of limits: %.2fx Sidereal (min=%.2f, max=%.2f)", 
		      absrate, MIN_RATE, MAX_RATE);
  } 
  //if (MountCode != 0xF0) {
    if (absrate > SKYWATCHER_LOWSPEED_RATE) {
      absrate = absrate / RAHighspeedRatio;
      useHighspeed = true;
    }
    //}
  period=(long)(((SKYWATCHER_STELLAR_DAY * (double)RAStepsWorm) / (double)RASteps360) / absrate);
  newstatus.direction = ((rate >= 0.0)? FORWARD: BACKWARD);
  //newstatus.slewmode=RAStatus.slewmode;
  newstatus.slewmode=SLEW;
  if (useHighspeed) newstatus.speedmode = HIGHSPEED; else newstatus.speedmode = LOWSPEED;
  if (RARunning) {
    if (newstatus.speedmode != RAStatus.speedmode)
      throw  EQModError(EQModError::ErrInvalidParameter, "Can not change rate while motor is running (speedmode differs).");
    if (newstatus.direction != RAStatus.direction)
      throw  EQModError(EQModError::ErrInvalidParameter, "Can not change rate while motor is running (direction differs).");
  }
  SetMotion(Axis1, newstatus);
  SetSpeed(Axis1, period);
}

void  Skywatcher::SetDERate(double rate)  throw (EQModError)
{
  double absrate=fabs(rate);
  unsigned long period=0;
  bool useHighspeed=false; 
  SkywatcherAxisStatus newstatus;

  if ((absrate < get_min_rate()) || (absrate > get_max_rate())) {
    throw  EQModError(EQModError::ErrInvalidParameter, "Speed rate out of limits: %.2fx Sidereal (min=%.2f, max=%.2f)", 
		      absrate, MIN_RATE, MAX_RATE);
  } 
  //if (MountCode != 0xF0) {
    if (absrate > SKYWATCHER_LOWSPEED_RATE) {
      absrate = absrate / DEHighspeedRatio;
      useHighspeed = true;
    }
    //}
  period=(long)(((SKYWATCHER_STELLAR_DAY * (double)DEStepsWorm) / (double)DESteps360) / absrate);
  newstatus.direction = ((rate >= 0.0)? FORWARD: BACKWARD);
  //newstatus.slewmode=DEStatus.slewmode;
  newstatus.slewmode=SLEW;
  if (useHighspeed) newstatus.speedmode = HIGHSPEED; else newstatus.speedmode = LOWSPEED;
  if (DERunning) {
    if (newstatus.speedmode != DEStatus.speedmode)
      throw  EQModError(EQModError::ErrInvalidParameter, "Can not change rate while motor is running (speedmode differs).");
    if (newstatus.direction != DEStatus.direction)
      throw  EQModError(EQModError::ErrInvalidParameter, "Can not change rate while motor is running (direction differs).");
  }
  SetMotion(Axis2, newstatus);
  SetSpeed(Axis2, period);
}

void Skywatcher::StartRATracking(double trackspeed) throw (EQModError)
{
  double rate;
  if (trackspeed != 0.0) rate = trackspeed / SKYWATCHER_STELLAR_SPEED;
  else rate = 0.0;
  if (rate != 0.0) {
    SetRARate(rate);
    if (!RARunning) StartMotor(Axis1);
  } else 
    StopMotor(Axis1);
}

void Skywatcher::StartDETracking(double trackspeed) throw (EQModError)
{
  double rate;
  if (trackspeed != 0.0) rate = trackspeed / SKYWATCHER_STELLAR_SPEED;
  else rate = 0.0;
  if (rate != 0.0) {
    SetDERate(rate);
    if (!DERunning) StartMotor(Axis2);
  } else 
    StopMotor(Axis2);
}

void Skywatcher::SetSpeed(SkywatcherAxis axis, unsigned long period) throw (EQModError)
{
  char cmd[7];
  long2Revu24str(period, cmd);
  if (axis==Axis1) RAPeriod=period; else DEPeriod=period;
  //IDLog("Setting axis %c speed to %d\n", axis, period);
  dispatch_command(SetStepPeriod, axis, cmd);
  read_eqmod(); 
}

void Skywatcher::SetTarget(SkywatcherAxis axis, unsigned long increment) throw (EQModError)
{
  char cmd[7];
  long2Revu24str(increment, cmd);
  //IDLog("Setting target for axis %c  to %d\n", axis, increment);
  dispatch_command(SetGotoTargetIncrement, axis, cmd);
  read_eqmod();  
}

void Skywatcher::SetTargetBreaks(SkywatcherAxis axis, unsigned long increment) throw (EQModError)
{
  char cmd[7];
  long2Revu24str(increment, cmd);
  //IDLog("Setting target for axis %c  to %d\n", axis, increment);
  dispatch_command(SetBreakPointIncrement, axis, cmd);
  read_eqmod();  
}


void Skywatcher::StartMotor(SkywatcherAxis axis) throw (EQModError)
{
  dispatch_command(StartMotion, axis, NULL);
  read_eqmod();
}

void Skywatcher::StopRA()  throw (EQModError)
{
  StopWaitMotor(Axis1);
}

void Skywatcher::StopDE()  throw (EQModError)
{
  StopWaitMotor(Axis2); 
}

void Skywatcher::SetMotion(SkywatcherAxis axis, SkywatcherAxisStatus newstatus) throw (EQModError)
{
  char motioncmd[3];
  SkywatcherAxisStatus *currentstatus;
 
  CheckMotorStatus(axis);
  if (axis == Axis1) currentstatus=&RAStatus; else currentstatus=&DEStatus;
  motioncmd[2]='\0';
  switch (newstatus.slewmode) {
  case SLEW: if (newstatus.speedmode == LOWSPEED) motioncmd[0] ='1'; else  motioncmd[0]='3'; break;
  case GOTO: if (newstatus.speedmode == LOWSPEED) motioncmd[0] ='2'; else  motioncmd[0]='0'; break;
  default: motioncmd[0] ='1'; break;
  }
  if (newstatus.direction == FORWARD) motioncmd[1] = '0'; else motioncmd[1] = '1';
  if ((newstatus.direction != currentstatus->direction) || (newstatus.speedmode != currentstatus->speedmode)
      || (newstatus.slewmode != currentstatus->slewmode))
    StopWaitMotor(axis);
  dispatch_command(SetMotionMode, axis, motioncmd);
  read_eqmod();
}


void Skywatcher::StopMotor(SkywatcherAxis axis)  throw (EQModError)
{
  dispatch_command(NotInstantAxisStop, axis, NULL);
  read_eqmod();
}

void Skywatcher::InstantStopMotor(SkywatcherAxis axis)  throw (EQModError)
{
  dispatch_command(InstantAxisStop, axis, NULL);
  read_eqmod();
}


void Skywatcher::StopWaitMotor(SkywatcherAxis axis)  throw (EQModError)
{
  bool *motorrunning;
  struct timespec wait;
  dispatch_command(NotInstantAxisStop, axis, NULL);
  read_eqmod();
  if (axis == Axis1) motorrunning=&RARunning; else motorrunning=&DERunning;
  wait.tv_sec=0; wait.tv_nsec=100000000; // 100ms
  ReadMotorStatus(axis);
  while (*motorrunning) {
    nanosleep(&wait, NULL);
    ReadMotorStatus(axis);
  }
}

/* Utilities */

void Skywatcher::CheckMotorStatus(SkywatcherAxis axis) throw (EQModError)
{
  struct timeval now;
 
  gettimeofday (&now, NULL);
  if (((now.tv_sec - lastreadmotorstatus[axis].tv_sec) + ((now.tv_usec - lastreadmotorstatus[axis].tv_usec)/1e6)) > SKYWATCHER_MAXREFRESH) 
    ReadMotorStatus(axis);
}

double Skywatcher::get_min_rate() {
  return MIN_RATE;
}

double Skywatcher::get_max_rate() {
  return MAX_RATE;
}

bool Skywatcher::dispatch_command(SkywatcherCommand cmd, SkywatcherAxis axis, char *command_arg) throw (EQModError)
{
  int err_code = 0, nbytes_written=0, nbytes_read=0;

  // Clear string
  command[0] = '\0';
  
  if (command_arg==NULL) 
    snprintf(command, SKYWATCHER_MAX_CMD, "%c%c%c%c", SkywatcherLeadingChar, cmd, axis, SkywatcherTrailingChar);
  else
    snprintf(command, SKYWATCHER_MAX_CMD, "%c%c%c%s%c", SkywatcherLeadingChar, cmd, axis, command_arg, SkywatcherTrailingChar);
  
  tcflush(fd, TCIOFLUSH);
  
  if  ( (err_code = tty_write_string(fd, command, &nbytes_written) != TTY_OK))
    {
      char ttyerrormsg[ERROR_MSG_LENGTH];
      tty_error_msg(err_code, ttyerrormsg, ERROR_MSG_LENGTH);
      throw EQModError(EQModError::ErrDisconnect, "tty write failed, check connection: %s", 
		       ttyerrormsg);
      return false;
   }

   return true;
}



bool Skywatcher::read_eqmod() throw (EQModError)
{
    int err_code = 0, nbytes_written=0, nbytes_read=0;

    // Clear string
    response[0] = '\0';

    //Have to onsider cases when we read ! (error) or 0x01 (buffer overflow)
    // Read until encountring a CR
    if ( (err_code = tty_read_section(fd, response, 0x0D, 15, &nbytes_read)) != TTY_OK)
    {
      char ttyerrormsg[ERROR_MSG_LENGTH];
      tty_error_msg(err_code, ttyerrormsg, ERROR_MSG_LENGTH);
      throw EQModError(EQModError::ErrDisconnect, "tty read failed, check connection: %s", 
		       ttyerrormsg);
      return false;
    }

    switch (response[0]) {
    case '=': break;
    case '!': throw  EQModError(EQModError::ErrCmdFailed,"Failed command %s-Response %s", command, response); return false;
    default:  throw  EQModError(EQModError::ErrInvalidCmd,"Invalid response to command %s-Response %s", command, response); return false;
    } 
    // Remove CR
    response[nbytes_read-1] = '\0';

    return true;
}

unsigned long Skywatcher::Revu24str2long(char *s) {
   unsigned long res = 0;
   res=HEX(s[4]); res <<= 4;
   res|=HEX(s[5]); res <<= 4;
   res|=HEX(s[2]); res <<= 4;
   res|=HEX(s[3]); res <<= 4;
   res|=HEX(s[0]); res <<= 4;
   res|=HEX(s[1]); 
   return res;
}


                
void Skywatcher::long2Revu24str(unsigned long n, char *str) {
char hexa[16] = {'0', '1', '2', '3', '4', '5', '6', '7', 
                 '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
   str[0]=hexa[(n & 0xF0) >> 4];
   str[1]=hexa[(n & 0x0F)];
   str[2]=hexa[(n & 0xF000) >> 12];
   str[3]=hexa[(n & 0x0F00) >> 8];
   str[4]=hexa[(n & 0xF00000) >> 20];
   str[5]=hexa[(n & 0x0F0000) >> 16];
   str[6]='\0';
}
