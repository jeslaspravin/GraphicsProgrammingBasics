// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

// For using Math constants
#define _USE_MATH_DEFINES 1

#include "targetver.h"

#include <cstdlib>
#include <stdio.h>
#include <tchar.h>

#include<fstream>
#include<iostream>
#include<assert.h>

#include "DataClasses/Vector.h"
#include "DataClasses/Matrix.h"

#include "Utils.h"

using namespace std;

#define GET_STREAM_DATA(msg,data) cout<<msg;cin>>data;
#define GET_TWO_STREAM_DATA(msg,data1,data2) cout<<msg;cin>>data1>>data2;


// TODO: reference additional headers your program requires here
