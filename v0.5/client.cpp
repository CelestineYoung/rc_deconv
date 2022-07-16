// ori netif.h
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/file.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>
#include <errno.h>
#include <signal.h>
#include <sys/wait.h>
// ori types.h
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/time.h>
// behave.cpp
#include <fstream>
// netif.cpp
#include <errno.h>
// test.cpp
#include <iomanip>

#define DEFAULT_PORT_NUMBER 6000
#define MAXMESG 2048
#define MAXARG 200

#define NULLCHAR '\000'

class Socket
{
public:
    int socketfd;
    struct sockaddr_in serv_addr;
};

Socket init_connection(char *host, int port);
int send_message(char *buf, Socket *sock);
int receive_message(char *buf, Socket *sock);
void close_connection(Socket *sock);

int wait_message(char *buf, Socket *sock);

typedef int Unum;  /* Uniform number           */
typedef int Pnum;  /* Position number          */
typedef int SPnum; /* Setplay position number  */

enum SenseType
{
    See_Msg,
    Hear_Msg,
    Sense_Msg
};

enum CMDType
{
    CMD_none,
    CMD_dash,
    CMD_turn,
    CMD_kick,
    CMD_catch,
    CMD_move,
    CMD_bye,
    CMD_change_view,
    CMD_turn_neck,
    CMD_say,
    CMD_sense_body
};

enum ObjType
{
    OBJ_Line,
    OBJ_Ball,
    OBJ_Marker,
    OBJ_Marker_Behind, /* Not seen */
    OBJ_Player
};

enum Vqual
{
    VQ_Low,
    VQ_High
};

enum Vwidth
{
    VW_Narrow,
    VW_Normal,
    VW_Wide
};

enum Kmode
{
    KO_Mine,
    KO_Theirs
};

enum Pmode
{
    PM_No_Mode,
    PM_Before_Kick_Off,
    PM_My_Kick_Off,
    PM_Their_Kick_Off,
    PM_My_Kick_In,
    PM_Their_Kick_In,
    PM_My_Corner_Kick,
    PM_Their_Corner_Kick,
    PM_My_Goal_Kick,
    PM_Their_Goal_Kick,
    PM_My_Free_Kick,
    PM_Their_Free_Kick,
    PM_My_Goalie_Free_Kick,    /* not a real play mode */
    PM_Their_Goalie_Free_Kick, /* not a real play mode */
    PM_Drop_Ball,
    PM_My_Offside_Kick,
    PM_Their_Offside_Kick,
    PM_Play_On,
    PM_Half_Time,
    PM_Time_Up,
    PM_Extended_Time
};

enum Bool
{
    FALSE,
    TRUE
};

enum SideLine
{
    SL_Left,
    SL_Right,
    SL_Top,
    SL_Bottom,

    SL_No_Line
};

enum MarkerType
{
    Goal_L,
    Goal_R,

    Flag_C,
    Flag_CT,
    Flag_CB,
    Flag_LT,
    Flag_LB,
    Flag_RT,
    Flag_RB,

    Flag_PLT,
    Flag_PLC,
    Flag_PLB,
    Flag_PRT,
    Flag_PRC,
    Flag_PRB,

    Flag_GLT,
    Flag_GLB,
    Flag_GRT,
    Flag_GRB,

    Flag_TL50,
    Flag_TL40,
    Flag_TL30,
    Flag_TL20,
    Flag_TL10,
    Flag_T0,
    Flag_TR10,
    Flag_TR20,
    Flag_TR30,
    Flag_TR40,
    Flag_TR50,

    Flag_BL50,
    Flag_BL40,
    Flag_BL30,
    Flag_BL20,
    Flag_BL10,
    Flag_B0,
    Flag_BR10,
    Flag_BR20,
    Flag_BR30,
    Flag_BR40,
    Flag_BR50,

    Flag_LT30,
    Flag_LT20,
    Flag_LT10,
    Flag_L0,
    Flag_LB10,
    Flag_LB20,
    Flag_LB30,

    Flag_RT30,
    Flag_RT20,
    Flag_RT10,
    Flag_R0,
    Flag_RB10,
    Flag_RB20,
    Flag_RB30,

    No_Marker
};

enum Ptype
{
    PT_None,
    PT_Goaltender,
    PT_Sweeper,
    PT_Defender,
    PT_Midfielder,
    PT_Forward
};

enum Pside
{
    PS_None,
    PS_Left,
    PS_Center,
    PS_Right
};

enum Fside
{ /* Side of the field */
  FS_Right,
  FS_Left
};

enum Utype
{
    UT_Defense,
    UT_Midfield,
    UT_Forward,
    UT_Left,
    UT_Center,
    UT_Right,
    UT_None
};

enum Ftype
{
    FT_None,
    FT_433,
    FT_442,
    FT_352,
    FT_72,
    FT_334,
    FT_244,
    FT_532,
    FT_right,
    FT_left
};

enum MCtype
{ /* Mark change type */
  MC_Obey,
  MC_Closest,
  MC_Open
};

enum HCtype
{ /* Home change type */
  HC_Obey,
  HC_Get_Open,
  HC_Shift,
  HC_Mark
};

enum SPAtype
{ /* Set play action type */
  SPA_None,
  SPA_Starter,
  SPA_Passer,  /* goes and passes   */
  SPA_Shooter, /* goes and shoots   */
  SPA_Knocker, /* aims at a point   */
  SPA_Blaster, /* blasts at a point */
  SPA_Getter   /* goes and ends setplay */
};

/* Communication targets and messages */

enum TargetType
{
    TT_Player,
    TT_Position,
    TT_Unit,
    TT_All
};

enum MessageType
{
    CMsg_new_coach, /* coach  */

    PMsg_none, /* player */
    PMsg_ping,
    PMsg_ping_ball,
    PMsg_ping_teammate,
    PMsg_ping_opponent,
    PMsg_ready_to_pass,
    PMsg_ready_to_receive,
    PMsg_passing_decision,
    PMsg_my_ball,
    PMsg_leaving_position,
    PMsg_already_there,
    PMsg_tired,
    PMsg_Setplay_Ready,
    PMsg_Setplay_OK_Ready,
    PMsg_Setplay_Starter,
    PMsg_Setplay_Ping_Starter,
    PMsg_marking,

    UMsg_assign_mark, /* unit   */
    UMsg_assign_position
};

/* Action Modes */

enum ActionMode
{
    AM_Unknown,

    AM_goaltend,

    AM_Localize,
    AM_Face_Ball,
    AM_Watch_Pass,
    AM_Recover,
    AM_Before_Kick_Off,
    AM_Setplay,
    AM_GetOnSide,

    AM_With_Ball,
    AM_Offense_Active,
    AM_Offense_Auxiliary,
    AM_Offense_Passive,

    AM_Defense_Active,
    AM_Defense_Auxiliary,
    AM_Defense_Passive
};

enum PassChoiceType
{
    PC_None,
    PC_Fixed,
    PC_Random,
    PC_DT_Max,
    PC_DT_Thresh,
    PC_Congestion
};

enum PassFilterType
{
    PF_None,
    PF_GoalDist_Congestion_And_Shot,
    PF_GoalDist_And_Congestion,
    PF_GoalDist_Or_Congestion,
    PF_XPos_Congestion_And_Shot,
    PF_XPos_And_Congestion,
    PF_XPos_Or_Congestion,
    PF_Breakaway,
    PF_BetterShot,
    PF_No_Opponent_Near
};

enum DodgeType
{
    DT_none,
    DT_all,
    DT_unless_with_ball,
    DT_only_with_ball
};

/* these are for things in kick.* */
typedef enum TURNDIR
{
    TURN_NONE = 0,
    TURN_CW = -1,
    TURN_CCW = 1,
    TURN_CLOSEST = 10,
    TURN_AVOID = 11 /* avoid any opponents */
} TurnDir;

typedef enum KICKTORES
{
    KT_None,
    KT_Success,
    KT_DidKick,
    KT_DidNothing,
    KT_TurnedToBall,
    KT_LostBall
} KickToRes;

typedef enum KICKMODE
{
    KM_None,
    KM_HardestKick,
    KM_Hard,
    KM_Moderate,
    KM_Quickly,
    KM_QuickestRelease
} KickMode;

/* these are for things in intercept.* */
typedef enum INTERCEPTRES
{
    BI_None,       /* no value yet */
    BI_Invalid,    /* could not get an answer */
    BI_Failure,    /* won;t be able to intercept ball */
    BI_CanChase,   /* we're getting there - returned a GoToPoint command*/
    BI_ReadyToKick /* ball is in kickable area, we haven;t done anything yet */
} InterceptRes;

typedef enum ACTIONQUEUERES
{
    AQ_ActionQueued,
    AQ_ActionNotQueued
} ActionQueueRes;

/* -*- Mode: C++ -*- */

/* utils.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

#define my_stamp printf("%d:%d.%d ", Mem->MyNumber, Mem->CurrentTime.t, Mem->CurrentTime.s);

int dump_core(char *);

void my_error(char *);
void my_error(char *, int);
void my_error(char *, int, int);
void my_error(char *, int, int, int);
void my_error(char *, int, int, int, int);
void my_error(char *, int, int, int, int, int);
void my_error(char *, int, int, int, int, int, int);
void my_error(char *, int, int, int, int, int, char, int);
void my_error(char *, float);
void my_error(char *, float, float);
void my_error(char *, float, float, float);
void my_error(char *, float, float, float, float);
void my_error(char *, float, int);
void my_error(char *, char *);
void my_error(char *, char, int, int);
void my_error(char *, char, int, float, float);

int closest_int(float x);

typedef float Value;
typedef float AngleRad;
typedef float AngleDeg;

inline AngleDeg Rad2Deg(AngleRad x) { return x * 180 / M_PI; }
inline AngleRad Deg2Rad(AngleDeg x) { return x * M_PI / 180; }

/* needed? */
/* inline float cos(AngleRad x) { return cos((float) x); } */
/* inline float sin(AngleRad x) { return sin((float) x); } */
/* inline float tan(AngleRad x) { return tan((float) x); } */

inline float Cos(AngleDeg x) { return cos(Deg2Rad(x)); }
inline float Sin(AngleDeg x) { return sin(Deg2Rad(x)); }
inline float Tan(AngleDeg x) { return tan(Deg2Rad(x)); }
inline AngleDeg ACos(float x) { return ((x) >= 1 ? 0 : ((x) <= -1 ? 180 : (Rad2Deg(acos(x))))); }
inline AngleDeg ASin(float x) { return ((x) >= 1 ? 90 : ((x) <= -1 ? -90 : (Rad2Deg(asin(x))))); }
inline AngleDeg ATan(float x) { return (Rad2Deg(atan(x))); }
inline AngleDeg ATan2(float x, float y) { return ((x == 0 && y == 0) ? 0 : (Rad2Deg(atan2(x, y)))); }

void NormalizeAngleDeg(int *);
void NormalizeAngleDeg(AngleDeg *);
void NormalizeAngleRad(AngleRad *);
AngleDeg GetNormalizeAngleDeg(AngleDeg);
float GetDistance(float *x, float *y, float *a, float *b);

#define FLOAT_EPS .001

#define Mod(a, b) (a - (b) * (int)((a) / (b)))
#define Sign(x) ((x) >= 0 ? 1 : -1)

float int_pow(float x, int p);
inline int Sqr(int x) { return x * x; }
inline float Sqr(float x) { return x * x; }
inline float Exp(float x, int y)
{
    float a = 1;
    for (int i = 0; i < y; i++)
        a *= x;
    return a;
}

inline float SumInfGeomSeries(float first_term, float r)
{
    return first_term / (1 - r);
}
float SumGeomSeries(float first_term, float r, int n);
/* returns -1 on error */
float SolveForLengthGeomSeries(float first_term, float r, float sum);
float SolveForFirstTermGeomSeries(float r, int n, float sum);
inline float SolveForFirstTermInfGeomSeries(float r, float sum)
{
    return sum * (1 - r);
}

#define signf(x) (((x) > 0.0) ? 1.0 : -1.0)
inline float Round(float x, int p = 0)
{
    x *= int_pow(10.0, -p);
    if (fmod(x, 1.0) >= .5)
        return ceil(x) / int_pow(10.0, -p);
    else
        return floor(x) / int_pow(10.0, -p);
}

const char char_for_num_array[16] =
    {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};

inline char char_for_num(int num)
{
    return char_for_num_array[num];
}

/* returns a pointer to a static buffer, so be careful! */
char *repeat_char(char c, int n);

class Time
{
public:
    int t; /* time from the server */
    int s; /* stopped clock cycles */

    Time(int vt = 0, int vs = 0)
    {
        t = vt;
        s = vs;
    }
    Time operator-(const int &a);
    int operator-(const Time &a);
    Time operator+(const int &a);
    int operator%(const int &a) { return (t + s) % a; }
    void operator-=(const int &a) { *this = *this - a; }
    void operator-=(const Time &a) { *this = *this - a; }
    void operator+=(const int &a) { *this = *this + a; }
    void operator++() { *this += 1; }
    void operator--() { *this -= 1; }
    Time operator=(const int &a)
    {
        t = a;
        s = 0;
        return *this;
    }
    bool operator==(const Time &a) { return (s == a.s) && (t == a.t); }
    bool operator==(const int &a) { return t == a; }
    bool operator!=(const Time &a) { return (s != a.s) || (t != a.t); }
    bool operator!=(const int &a) { return t != a; }
    bool operator<(const Time &a) { return (t < a.t) || (t == a.t && s < a.s); }
    bool operator<(const int &a) { return t < a; }
    bool operator<=(const Time &a) { return (t < a.t) || (t == a.t && s <= a.s); }
    bool operator<=(const int &a) { return t <= a; }
    bool operator>(const Time &a) { return (t > a.t) || (t == a.t && s > a.s); }
    bool operator>(const int &a) { return t > a; }
    bool operator>=(const Time &a) { return (t > a.t) || (t == a.t && s >= a.s); }
    bool operator>=(const int &a) { return t >= a; }
    bool operator!() { return (s == 0) && (t == 0); }

    Bool CanISubtract(const Time &a);
};

#define Min(x, y) ((x) < (y) ? (x) : (y))
#define Max(x, y) ((x) > (y) ? (x) : (y))
#define MinMax(min, x, max) Min(Max((min), (x)), (max))

double get_double(char **str_ptr);
double get_double(char *str);
float get_float(char **str_ptr);
float get_float(char *str);
int get_int(char **str_ptr);
int get_int(char *str);
void get_word(char **str_ptr);
void get_next_word(char **str_ptr);
void get_token(char **str_ptr);
void advance_to(char c, char **str_ptr);
void advance_past_space(char **str_ptr);

int put_int(char *str, int num);
int put_float(char *str, float fnum, int precision);

void BubbleSort(int length, int *elements, float *keys);
int BinarySearch(int length, float *elements, float key);
void StrReplace(char *str, char oldchar, char newchar);

int int_random(int n);
float range_random(float lo, float hi);
int very_random_int(int n);

float weighted_avg(float val1, float val2, float w1, float w2);

void GetStampedName(char *name, char *outputName);

/* -*- Mode: C++ -*- */

/* geometry.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        Vector Class                                          */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

class Vector
{
public:
    float x;
    float y;

    Vector(float vx = 0, float vy = 0)
    {
        x = vx;
        y = vy;
    }
    Vector operator-() { return Vector(-x, -y); }
    Vector operator+(const Vector &a) { return Vector(x + a.x, y + a.y); }
    Vector operator-(const Vector &a) { return Vector(x - a.x, y - a.y); }
    Vector operator*(const float &a) { return Vector(x * a, y * a); }
    Vector operator*(const Vector &a) { return Vector(x * a.x, y * a.y); }
    Vector operator/(const float &a) { return Vector(x / a, y / a); }
    Vector operator/(const Vector &a) { return Vector(x / a.x, y / a.y); }
    void operator=(const float &a)
    {
        x = a;
        y = a;
    }
    void operator+=(const Vector &a)
    {
        x += a.x;
        y += a.y;
    }
    void operator+=(const float &a)
    {
        x += a;
        y += a;
    }
    void operator-=(const Vector &a)
    {
        x -= a.x;
        y -= a.y;
    }
    void operator-=(const float &a)
    {
        x -= a;
        y -= a;
    }
    void operator*=(const float &a)
    {
        x *= a;
        y *= a;
    }
    void operator/=(const float &a)
    {
        x /= a;
        y /= a;
    }
    bool operator!=(const Vector &a) { return (x != a.x) || (y != a.y); }
    bool operator!=(const float &a) { return (x != a) || (y != a); }
    bool operator==(const Vector &a) { return (x == a.x) && (y == a.y); }
    friend std::ostream &operator<<(std::ostream &os, Vector v) { return os << "(" << v.x << ", " << v.y << ")"; }
    void Print() { std::cout << *this << std::endl; }
    void ins(float vx = 0, float vy = 0)
    {
        x = vx;
        y = vy;
    }
    float mod() { return sqrt(x * x + y * y); }
    float mod2() { return x * x + y * y; }
    float dist(const Vector &a) { return (*this - a).mod(); }
    float dist2(const Vector &a) { return (*this - a).mod2(); }
    Vector Normalize() { return SetLength(1.0); }
    Vector SetLength(float len) { return (*this) * (len / mod()); }
    AngleDeg dir() { return ATan2(y, x); }
    Vector rotate(AngleDeg angle)
    {
        return Vector(this->mod() * Cos(this->dir() + angle), this->mod() * Sin(this->dir() + angle));
    }
    Vector Global2Relative(Vector orig, float ang);
    Vector Relative2Global(Vector orig, float ang);
    Bool InFrontOf(const Vector &a) { return x > a.x ? TRUE : FALSE; }
    Bool InFrontOf(const float &a) { return x > a ? TRUE : FALSE; }
    Bool Behind(const Vector &a) { return x < a.x ? TRUE : FALSE; }
    Bool Behind(const float &a) { return x < a ? TRUE : FALSE; }

    Bool ApproxEqual(const Vector &a)
    {
        return (fabs(x - a.x) < FLOAT_EPS && fabs(y - a.y) < FLOAT_EPS) ? TRUE : FALSE;
    }

    Vector WeightedAverage(const Vector &a, float w)
    {
        return Vector((1.0 - w) * x + w * a.x, (1.0 - w) * y + w * a.y);
    }
};

inline Vector Polar2Vector(float mod, AngleDeg ang)
{
    return Vector(mod * Cos(ang), mod * Sin(ang));
}

inline Vector Vector::Global2Relative(Vector orig, float ang)
{
    return (*this - orig).rotate(-ang);
}

inline Vector Vector::Relative2Global(Vector orig, float ang)
{
    return (rotate(ang) + orig);
}

inline Vector operator*(const float &a, const Vector &v)
{
    return Vector(v.x * a, v.y * a);
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                       Ray Class                                              */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/* stores origin and direction (which is normalized) */
class Line;
class Rectangle;
class Ray
{
public:
    Ray()
    {
        origin = Vector(0, 0);
        direction = Vector(1, 0);
    }
    Ray(Vector orig, Vector dir);
    Ray(Vector orig, float ang)
    {
        origin = orig;
        direction = Polar2Vector(1.0, ang);
    }

    Bool OnRay(Vector pt);
    Bool InRightDir(Vector pt); // more lenient than above about distance off ray

    Bool intersection(Ray r, Vector *pPt);
    Bool intersection(Line l, Vector *pPt);

    int CircleIntersect(float rad, Vector center,
                        Vector *psol1, Vector *psol2);

    Vector GetClosestPoint(Vector pt);

    Vector RectangleIntersection(Rectangle R);

protected:
    friend Line;
    friend Rectangle;
    Vector origin;
    Vector direction;
};

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                      Line Class                                              */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
class Line
{
public:
    Line(float x_coef = 0.0, float y_coef = 0.0, float constant = 0.0);
    Line(Ray r)
    {
        LineFromRay(r);
    }

    void LineFromTwoPoints(Vector pt1, Vector pt2);
    void LineFromRay(Ray r);
    void LineFromRay(Vector orig, Vector dir)
    {
        LineFromRay(Ray(orig, dir));
    }
    void LineFromRay(Vector orig, AngleDeg dir)
    {
        LineFromRay(orig, Polar2Vector(1.0, dir));
    }

    Bool PointOnLine(float x, float y);
    inline Bool PointOnLine(Vector pt) { return PointOnLine(pt.x, pt.y); }

    float dist(Vector pt);
    float dist2(Vector pt);
    float angle(); /* returns the angle of the line in [-90, 90] */

    Bool InBetween(Vector pt, Vector end1, Vector end2);
    Vector GetClosestPtInBetween(Vector pt, Vector end1, Vector end2);

    /* the buffer should really be linked to an option or something,
       but 1.0 is okay */
    inline Bool OnLine(Vector pt, float buffer = 1.0)
    {
        return (dist2(pt) < Sqr(buffer)) ? TRUE : FALSE;
    }

    Vector ProjectPointUsingCircle(Vector pt);
    inline Vector ProjectPoint(Vector pt) { return intersection(perpendicular(pt)); }

    float get_y(float x);
    float get_x(float y);
    Vector intersection(Line l);
    inline Line perpendicular(Vector pt) { return Line(B, -A, A * pt.y - B * pt.x); }

    Bool RayIntersection(Ray r, Vector *ppt);

    Line shift_y(float val)
    {
        return Line(A, B, C - val * B);
    }
    Line shift_x(float val)
    {
        return Line(A, B, C - val * A);
    }

    /* returns whether the projection of pt1 is closer to targ_pt than the
       projection of pt2 */
    Bool IsPtCloserToPtOnLine(Vector pt1, Vector pt2, Vector targ_pt);

    Bool HalfPlaneTest(Vector pt); // return TRUE on top/left part of plane
    Bool SameSlope(Line l);

    friend std::ostream &operator<<(std::ostream &os, Line l)
    {
        return os << "#L(" << l.A << ", " << l.B << ", " << l.C << ")";
    }
    void Print(void) { std::cout << *this << std::endl; }

    // private:
    float A, B, C; /* the three coeffs in the line equation */
                   /* Ax + By + C = 0 */
};

inline Line LineFromTwoPoints(Vector pt1, Vector pt2)
{
    Line l;
    l.LineFromTwoPoints(pt1, pt2);
    return l;
}
inline Line LineFromRay(Ray r)
{
    Line l;
    l.LineFromRay(r);
    return l;
}
inline Line LineFromRay(Vector orig, Vector dir)
{
    Line l;
    l.LineFromRay(orig, dir);
    return l;
}
inline Line LineFromRay(Vector orig, AngleDeg dir)
{
    Line l;
    l.LineFromRay(orig, dir);
    return l;
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        Rectangle Class                                       */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

class Rectangle
{
public:
    Rectangle();
    Rectangle(const float left, const float right,
              const float top, const float bottom);
    Rectangle(const Vector center, const Vector size);

    bool operator==(const Rectangle &a)
    {
        return (left_x == a.left_x) && (right_x == a.right_x) &&
               (top_y == a.top_y) && (bottom_y == a.bottom_y);
    }

    inline float TopY() { return top_y; }
    inline float BottomY() { return bottom_y; }
    inline float RightX() { return right_x; }
    inline float LeftX() { return left_x; }

    inline Vector TopLeftCorner() { return Vector(left_x, top_y); }
    inline Vector TopRightCorner() { return Vector(right_x, top_y); }
    inline Vector BottomLeftCorner() { return Vector(left_x, bottom_y); }
    inline Vector BottomRightCorner() { return Vector(right_x, bottom_y); }

    inline Line LeftEdge() { return LineFromTwoPoints(TopLeftCorner(), BottomLeftCorner()); }
    inline Line RightEdge() { return LineFromTwoPoints(TopRightCorner(), BottomRightCorner()); }
    inline Line TopEdge() { return LineFromTwoPoints(TopLeftCorner(), TopRightCorner()); }
    inline Line BottomEdge() { return LineFromTwoPoints(BottomLeftCorner(), BottomRightCorner()); }

    Line GetEdge(int n);    // order: top, right, bot, left; starting at 0
    Vector GetPoint(int n); // order: TL, TR, BR, BL; starting at 0

    inline float Width() { return right_x - left_x; }
    inline float Height() { return bottom_y - top_y; }

    inline Vector Center() { return Vector(LeftX() + Width() / 2, TopY() + Height() / 2); }

    Bool IsWithin(const Vector &p);
    Vector nearestHEdge(const Vector &p);
    Vector nearestVEdge(const Vector &p);
    Vector nearestEdge(const Vector &p);

    float DistanceToEdge(const Vector &p);
    Vector AdjustToWithin(const Vector &p);

    Line nearestHEdgeLine(const Vector &p);
    Line nearestVEdgeLine(const Vector &p);
    Line nearestEdgeLine(const Vector &p);

    Vector random();

    Rectangle shrink(float val)
    {
        return expand(-val);
    }
    Rectangle expand(float val);

    Rectangle expandLeft(float val)
    {
        return Rectangle(left_x - val, right_x, top_y, bottom_y);
    }
    Rectangle expandRight(float val)
    {
        return Rectangle(left_x, right_x + val, top_y, bottom_y);
    }
    Rectangle expandTop(float val)
    {
        return Rectangle(left_x, right_x, top_y - val, bottom_y);
    }
    Rectangle expandBottom(float val)
    {
        return Rectangle(left_x, right_x, top_y, bottom_y + val);
    }
    Rectangle shrinkLeft(float val)
    {
        return expandLeft(-val);
    }
    Rectangle shrinkRight(float val)
    {
        return expandRight(-val);
    }
    Rectangle shrinkTop(float val)
    {
        return expandTop(-val);
    }
    Rectangle shrinkBottom(float val)
    {
        return expandBottom(-val);
    }

    friend std::ostream &operator<<(std::ostream &os, Rectangle r)
    {
        return os << "RECTANGLE:  x = " << r.LeftX() << " to " << r.RightX()
                  << "   y = " << r.TopY() << " to " << r.BottomY();
    }
    void Print();

    Vector RayIntersection(Ray r);

private:
    float left_x;
    float right_x;
    float top_y;
    float bottom_y;
};

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                  Miscellaneous Geometry Functions                            */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

int QuadraticFormula(float a, float b, float c, float *psol1, float *psol2);

inline int RayCircleIntersect(Ray r, float rad, Vector center,
                              Vector *psol1, Vector *psol2)
{
    return r.CircleIntersect(rad, center, psol1, psol2);
}

int LineCircleIntersect(Line l, float rad, Vector center,
                        Vector *psol1, Vector *psol2);

Vector AdjustPtToRectOnLine(Vector pt, Rectangle r, Line l);

Bool InBetween(Vector pt, Vector end1, Vector end2);

Vector PointInBetween(Vector pt1, Vector pt2, float pt1dist = 0);

AngleDeg AngleBisect(AngleDeg a1, AngleDeg a2);

Vector GetClosestPtInBetween(Vector pt, Vector end1, Vector end2);

Bool IsPointInCone(Vector pt, float wid_dist_ratio, Vector end, Vector vert);

/* -*- Mode: C++ -*- */
/* MemOption.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

struct option_t
{
    char optname[50];
    void *vptr;
    int vsize;
};

enum InputType
{
    V_INT,
    V_FLOAT,
    V_BOOL,
    V_STRING,
    V_ONOFF,
    V_NONE
};

#define MAX_TEAMNAME_LEN 20
#define MAX_HOST_LEN 50
#define MAX_FILE_LEN 256
#define MAX_FP_LEN 20
#define MAX_TREE_STEM_LEN 50

/* Things to be read at startup that never change */
class OptionInfo
{
public:
    OptionInfo();
    void GetOptions(int, char **);
    void GetOptions(char *);

    /* Client version params */
    Bool VP_test_l;
    Bool VP_test_r;
    Bool VP_test;
    Bool VP_train_DT;
    Bool VP_use_DT;

    /* Initialize params */
    int IP_my_score;
    int IP_their_score;
    int IP_reconnect; /* If non-zero, reconnect to that unum */
    /* Client params */
    char MyTeamName[MAX_TEAMNAME_LEN];
    Bool CP_goalie;
    Bool CP_save_log;
    int CP_save_freq;
    Bool CP_save_sound_log;
    int CP_save_sound_freq;
    int CP_save_action_log_level;
    int CP_save_action_freq;
    float CP_send_ban_recv_step_factor;
    int CP_interrupts_per_cycle;
    int CP_interrupts_left_to_act;
    float CP_max_conf;
    float CP_min_valid_conf;
    float CP_conf_decay;
    float CP_player_conf_decay;
    float CP_ball_conf_decay;
    float CP_max_player_move_factor; /* multiply by SP_player_speed_max to see how far a player can be
                        from its last position and still be considered the same player */
    int CP_max_say_interval;
    float CP_ball_moving_threshold;
    float CP_dodge_angle_buffer;
    float CP_dodge_distance_buffer;
    float CP_dodge_power;
    float CP_dodge_angle;
    char CP_tree_stem[MAX_TREE_STEM_LEN];
    int CP_DT_evaluate_interval;
    int CP_say_tired_interval;
    float CP_tired_buffer;
    Bool CP_set_plays;
    int CP_Setplay_Delay;
    int CP_Setplay_Say_Delay;
    int CP_Setplay_Max_Delay;
    int CP_Setplay_Time_Limit;
    float CP_kickable_buffer;
    int CP_mark_persist_time;
    float CP_track_max_distance;
    float CP_track_min_distance;
    Bool CP_pull_offsides;
    Bool CP_pull_offsides_when_winning;
    Bool CP_spar;
    Bool CP_mark;
    Bool CP_communicate;
    int CP_change_view_for_ball_cycles;
    float CP_defer_kick_to_teammate_buffer;
    float CP_scan_overlap_angle;

    float CP_pull_offside_threshold;
    float CP_pull_offside_buffer;

    float CP_ball_forget_angle_buf;
    float CP_player_forget_angle_buf;
    float CP_ball_forget_dist_buf;
    float CP_player_forget_dist_buf;

    /* pat added these */
    Bool CP_use_new_position_based_vel;
    Bool CP_stop_on_error;

    /* these parameters affect turnball */
    float CP_max_turn_kick_pow;
    float CP_opt_ctrl_dist;
    float CP_closest_margin;
    float CP_dokick_factor;

    /* these basically affect the way turnball stops */
    float CP_KickTo_err;
    float CP_max_ignore_vel;

    int CP_kick_time_space;
    float CP_max_est_err;
    float CP_holdball_kickable_buffer;
    int CP_stop_ball_power;
    int CP_possessor_intercept_space;
    int CP_can_keep_ball_cycle_buffer;

    /* no longer used
    float CP_hard_kick_margin;
    float CP_hard_kick_factor;
    float CP_hard_kick_end_turn_dist; */
    float CP_hard_kick_dist_buffer;
    int CP_max_hard_kick_angle_err;
    /* angle off perpendicualr to start ball for hardest kick */
    int CP_hardest_kick_ball_ang;
    float CP_hardest_kick_ball_dist;
    int CP_hardest_kick_player_ang;
    float CP_max_dash_help_kick_angle;

    int CP_max_go_to_point_angle_err;
    int CP_max_int_lookahead;
    float CP_intercept_close_dist;
    int CP_intercept_step;
    int CP_my_intercept_step;
    int CP_intercept_aim_ahead;
    int CP_no_turn_max_cyc_diff;    /* used for normal interception */
    float CP_no_turn_max_dist_diff; /* used for ball_path intercept */
    float CP_turnball_opp_worry_dist;
    float CP_collision_buffer;
    float CP_behind_angle;
    int CP_time_for_full_rotation;
    float CP_ball_vel_invalidation_factor;

    /* dribble params */
    int CP_dribble_dash_pow;
    float CP_dribble_ball_dist;
    /* dist where opponent starts to affect where we dribble ball */
    float CP_dribble_ignore_opp_dist;
    /* dist of opponent that makes us go to DM_Strict mode */
    float CP_dribble_worry_opp_dist;
    /* angle we normally like to dribble at */
    float CP_dribble_angle_norm;
    /* max and min distnaces to worry about dodging a player */
    float CP_dribble_dodge_max_dist;
    /* angle diff to make us turn if dodging */
    float CP_dribble_dodge_angle_err;
    /* how far off in expected angle we let a ball before we kick it to correct */
    float CP_dribble_exp_angle_buffer;
    /* if drib_ang > 180 - X, we will just dribble on the side where the ball is */
    float CP_dribble_angle_ignore_buffer;
    float CP_dribble_dodge_close_dist;
    float CP_can_dribble_cone_ratio;
    float CP_dribble_towards_length;
    float CP_dribble_sideline_buffer;
    float CP_dribble_circle_inner_rad;
    float CP_dribble_circle_outer_rad;
    float CP_dribble_circle_ang; // angle realtive to dribble angle to look for players
    Bool CP_dribble_scan_field;

    float CP_move_imp_1v1_initial;
    float CP_move_imp_1v1_inc;
    float CP_move_imp_1v1_threshold;
    float CP_at_point_buffer;
    float CP_overrun_dist;
    float CP_def_block_dist;
    float CP_def_block_dist_ratio;
    float CP_overrun_buffer;
    float CP_breakaway_buffer;
    float CP_our_breakaway_kickable_buffer;
    float CP_their_breakaway_front_kickable_buffer;
    float CP_their_breakaway_back_kickable_buffer;
    float CP_goalie_breakaway_kickable_buffer;

    float CP_breakaway_approach_x;
    float CP_breakaway_approach_y;
    int CP_breakaway_targ_valid_time;
    int CP_breakaway_min_goalie_steal_time;
    int CP_breakaway_kick_run_min_cycles;
    int CP_breakaway_kick_run_max_cycles;
    float CP_their_breakaway_min_cone_dist_wid;
    float CP_our_breakaway_min_cone_dist_wid;
    float CP_breakaway_middle_buffer;
    float CP_breakaway_kick_run_worry_dist;
    int CP_breakaway_mode; // used to test diff breakaway styles

    float CP_beat_offsides_buffer;
    float CP_beat_offsides_threshold;
    float CP_beat_offsides_max_x;
    float CP_congestion_epsilon;
    float CP_back_pass_opponent_buffer;
    float CP_back_pass_offside_buffer;
    float CP_min_less_congested_pass_dist;

    float CP_cycles_to_kick;

    /* parameters for moving to a standing ball */
    float CP_static_kick_dist_err;
    float CP_static_kick_ang_err;
    /* no longer used
    float CP_static_kick_dist;
    float CP_static_kick_ang;
    float CP_static_kick_overrun_dist;
    */

    float CP_goalie_baseline_buffer;
    float CP_goalie_scan_angle_err;
    float CP_goalie_at_point_buffer;
    float CP_goalie_vis_angle_err;
    float CP_goalie_max_shot_distance;
    float CP_goalie_min_pos_dist;
    float CP_goalie_max_pos_dist;
    float CP_goalie_max_forward_percent;
    float CP_goalie_ball_ang_for_corner;
    float CP_goalie_max_come_out_dist;
    float CP_goalie_ball_dist_for_corner;
    float CP_goalie_ball_dist_for_center;
    float CP_goalie_free_kick_dist;
    float CP_goalie_go_to_ball_cone_ratio;
    int CP_goalie_warn_space;
    Bool CP_goalie_comes_out;
    int CP_goalie_catch_wait_time;
    float CP_goalie_opponent_dist_to_block;
    float CP_goalie_position_weight_dist;
    int CP_goalie_narrow_sideline_cyc;
    float CP_goalie_no_buffer_dist;

    float CP_clear_ball_ang_step;
    float CP_clear_ball_cone_ratio;
    float CP_clear_ball_max_dist;
    float CP_clear_offensive_min_horiz_dist;
    float CP_clear_offensive_min_angle;

    float CP_should_cross_corner_dist;
    float CP_should_cross_baseline_buffer;
    float CP_should_move_to_cross_corner_dist;
    float CP_cross_pt_x;
    float CP_cross_pt_y;
    float CP_cross_target_vel;

    float CP_dont_dribble_to_middle_min_x;

    /* not used anymore
      float CP_hardest_kick_shot_distance;
      float CP_moderate_kick_shot_distance;
    */
    float CP_good_shot_distance;
    float CP_shot_distance;
    int CP_cycles_to_kick_buffer;
    float CP_shot_speed;
    int CP_shot_goalie_react_buffer;
    int CP_good_shot_goalie_react_buffer;
    int CP_better_shot_cyc_diff;
    // float CP_breakaway_shot_distance; no longer used

    /* Formation params */
    char FP_initial_formation[MAX_FP_LEN];
    char FP_formation_when_tied[MAX_FP_LEN];
    char FP_formation_when_losing[MAX_FP_LEN];
    char FP_formation_when_losing_lots[MAX_FP_LEN];
    char FP_formation_when_winning[MAX_FP_LEN];
    char FP_initial_hc_method[MAX_FP_LEN];
    char FP_initial_mc_method[MAX_FP_LEN];
    int FP_initial_player_1_pos;
    int FP_initial_player_2_pos;
    int FP_initial_player_3_pos;
    int FP_initial_player_4_pos;
    int FP_initial_player_5_pos;
    int FP_initial_player_6_pos;
    int FP_initial_player_7_pos;
    int FP_initial_player_8_pos;
    int FP_initial_player_9_pos;
    int FP_initial_player_10_pos;
    int FP_initial_player_11_pos;
    int FP_goalie_number;

    /* Server params */
    float SP_pitch_length;
    float SP_pitch_width;
    float SP_pitch_margin;
    float SP_penalty_area_length;
    float SP_penalty_area_width;
    float SP_goal_area_length;
    float SP_goal_area_width;
    float SP_penalty_spot_dist;
    float SP_corner_arc_r;
    float SP_free_kick_buffer;
    int SP_after_goal_wait;
    float SP_feel_distance;
    int SP_num_lines;
    int SP_num_markers;
    float SP_unum_far_length;
    float SP_unum_too_far_length;
    float SP_team_far_length;
    float SP_team_too_far_length;

    float SP_version;
    int SP_team_size;
    int SP_half;
    char SP_host[MAX_HOST_LEN];
    float SP_goal_width;
    float SP_player_size;
    float SP_player_decay;
    float SP_player_rand;
    float SP_player_weight;
    float SP_player_speed_max;
    float SP_stamina_max;
    float SP_stamina_inc;
    float SP_recover_dec_thr;
    float SP_recover_min;
    float SP_recover_dec;
    float SP_effort_dec_thr;
    float SP_effort_min;
    float SP_effort_dec;
    float SP_effort_inc_thr;
    float SP_effort_inc;
    float SP_ball_size;
    float SP_ball_decay;
    float SP_ball_rand;
    float SP_ball_weight;
    float SP_ball_speed_max;
    float SP_dash_power_rate;
    float SP_kick_power_rate;
    float SP_kickable_margin;
    float SP_kickable_area;
    float SP_catch_prob;
    float SP_catch_area_l;
    float SP_catch_area_w;
    float SP_max_power;
    float SP_min_power;
    float SP_max_moment;
    float SP_min_moment;
    float SP_max_neck_angle;
    float SP_min_neck_angle;
    float SP_max_neck_moment;
    float SP_min_neck_moment;
    float SP_visible_angle;
    float SP_visible_dist;
    float SP_audio_cut_dist;
    float SP_dist_qstep;
    float SP_land_qstep;
    float SP_ckmargin;
    float SP_wind_dir;
    float SP_wind_force;
    float SP_wind_rand;
    Bool SP_wind_none;
    Bool SP_wind_random;
    int SP_half_time;
    int SP_port;
    int SP_coach_port;
    int SP_olcoach_port;
    int SP_simulator_step;
    int SP_send_step;
    int SP_recv_step;
    int SP_say_msg_size;
    int SP_hear_max;
    int SP_hear_inc;
    int SP_hear_decay;
    int SP_catch_ban_cycle;
    Bool SP_coach_mode;
    Bool SP_coach_w_referee_mode;
    int SP_say_coach_cnt_max;
    int SP_say_coach_msg_size;
    int SP_send_vi_step;
    int SP_look_step;

    Bool SP_use_offside;
    Bool SP_forbid_kickoff_offside;
    char SP_logfile[MAX_FILE_LEN];
    char SP_recfile[MAX_FILE_LEN];
    Bool SP_rec_log;
    int SP_rec_ver;
    char SP_replay[MAX_FILE_LEN];
    Bool SP_verbose;
    Bool SP_send_log;
    float SP_offside_area;
    float SP_inertia_moment;
    int SP_sense_body_step;
    float SP_offside_kick_margin;
    Bool SP_record_messages;
};

/* -*- Mode: C++ -*- */

/* MemPlayer.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

class Command
{
public:
    CMDType type;

    float power;
    float angle;
    float x;
    float y;
    Vqual qual;
    Vwidth width;

    char command[MAXMESG];
    Time time;

    Command() { type = CMD_none; }

    inline Bool valid() { return (Bool)(type != CMD_none); }
    inline Bool valid(Time t) { return (Bool)(type != CMD_none && t == time); }
};

/*****************************************/
/*****************************************/
/*****************************************/

class PlayerInfo : public OptionInfo
{
public:
    PlayerInfo();
    virtual ~PlayerInfo();
    void Initialize();

    void SetPlayMode(Pmode mode);
    void sanitize_time(Time &tm);

    void EstimateMyPos();
    void EstimateMyVel(Time time);
    Vector NewVelFromDash(Vector old_vel, float dash_power);
    void UpdateFromMyAction(Time time);
    void update_self_estimate(Time time);
    void update_self_neck_rel_ang(Time time);
    void update_stamina(Time time);
    void reset_stamina();
    Time update_time(int time);

    virtual void VerifyDash(float *dash_power);

    /* You can specify a flag at compile time so that all the LogAction calls
       disappear
       If you add a LogAction, make sure to add it in both places */
    /* The # at the end of the function name is the total # of args */
#ifdef NO_ACTION_LOG
    inline void nothing_func(void)
    {
    }
#define LogAction1(x1) nothing_func()
#define LogAction2(x1, x2) nothing_func()
#define LogAction3(x1, x2, x3) nothing_func()
#define LogAction4(x1, x2, x3, x4) nothing_func()
#define LogAction5(x1, x2, x3, x4, x5) nothing_func()
#define LogAction6(x1, x2, x3, x4, x5, x6) nothing_func()
#define LogAction7(x1, x2, x3, x4, x5, x6, x7) nothing_func()
#define LogAction8(x1, x2, x3, x4, x5, x6, x7, x8) nothing_func()
#define LogAction9(x1, x2, x3, x4, x5, x6, x7, x8, x9) nothing_func()
#define LogAction10(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10) nothing_func()
/*
  inline void LogAction(int, char*) {}
  inline void LogAction(int, char*, char*) {}
  inline void LogAction(int, char*, char*, char*) {}
  inline void LogAction(int, char*, char) {}
  inline void LogAction(int, char*, float) {}
  inline void LogAction(int, char*, float, int) {}
  inline void LogAction(int, char*, float, float) {}
  inline void LogAction(int, char*, float, float, float) {}
  inline void LogAction(int, char*, float, float, float, float) {}
  inline void LogAction(int, char*, float, float, float, float, float) {}
  inline void LogAction(int, char*, float, float, float, float, float, float) {}
  inline void LogAction(int, char*, int) {}
  inline void LogAction(int, char*, int, int) {}
  */
#else
    void LogAction2(int level, char *str); /* will be timestamped automatically */
    void LogAction3(int level, char *str, char *param);
    void LogAction4(int level, char *str, char *param1, char *param2);
    void LogAction3(int level, char *str, char c1);
    void LogAction3(int level, char *str, float f1);
    void LogAction4(int level, char *str, float f1, int d1);
    void LogAction4(int level, char *str, float f1, float f2);
    void LogAction5(int level, char *str, float f1, float f2, float f3);
    void LogAction6(int level, char *str, float f1, float f2, float f3, float f4);
    void LogAction7(int level, char *str, float f1, float f2, float f3, float f4, float f5);
    void LogAction8(int level, char *str, float f1, float f2, float f3, float f4, float f5, float f6);
    void LogAction3(int level, char *str, int d1);
    void LogAction4(int level, char *str, int d1, int d2);
    void LogAction4(int level, char *str, int d1, float f1);
    void LogAction5(int level, char *str, int d1, float f1, float f2);
    void LogAction6(int level, char *str, int d1, float f1, float f2, float f3);
    void LogAction7(int level, char *str, int d1, float f1, float f2, float f3, float f4);

    void LogAction7(int level, char *str, int d1, int d2, float f1, float f2, float f3);
#endif

    Socket *sock;

    char SaveLogFileName[MAX_FILE_LEN];
    FILE *SaveLogFile;
    int SaveLogCounter;
    char SaveSoundLogFileName[MAX_FILE_LEN];
    FILE *SaveSoundLogFile;
    int SaveSoundLogCounter;
    char SaveActionLogFileName[MAX_FILE_LEN];
    FILE *SaveActionLogFile;
    int SaveActionLogCounter;

    char MySide;                          // 'l' or 'r'
    char TheirSide;                       // 'l' or 'r'
    char TheirTeamName[MAX_TEAMNAME_LEN]; // The name of their team
    int MyTeamNameLen;                    // strlen of MyTeamName
    Unum MyNumber;                        // uniform number
    int Initialized;
    Bool ServerAlive; // is the server going?
    Bool TestVersion; // is this the test version of the program?
    Bool CoachActive; // is there an on-line coach?

    int TimerInterval;
    Bool ClockStopped;
    int StoppedClockMSec; // time elapsed in before_kick_off_mode (msec)

    Time CurrentTime;   // Most recently observed time
    Time LastTime;      // previous time
    Time LastSightTime; // last time I saw
    Time LastSenseTime; // last time I sensed
    Time LastSoundTime; // last time I heard

    /* Don't need these after server 5.23
    void SetRealTimeOfLastSend(struct timeval& tv)
      { real_time_of_last_send = tv;}
    struct timeval GetRealTimeOfLastSend()
      { return real_time_of_last_send;}
    Bool TooSoonForAnotherSend();
    */

    int LastSightInterval;                                                        // cycles between last 2 sights
    inline Time PreviousSightTime() { return LastSightTime - LastSightInterval; } // previous time I saw

    Command *Action;
    Command *LastAction;
    Command ChangeView;
    Command TurnNeck;
    inline Bool ResendNeeded()
    {
        return (RequestResend && CurrentTime == LastActionTime() && ResendTime == CurrentTime && LastActionType() == ResendType) ? TRUE : FALSE;
    }
    Bool RequestResend;
    CMDType ResendType;
    Time ResendTime;

    inline AngleDeg LastActionAngle() { return LastAction->angle; }
    inline float LastActionPower() { return LastAction->power; }
    inline float LastActionX() { return LastAction->y; }
    inline float LastActionY() { return LastAction->x; }
    inline CMDType LastActionType() { return LastAction->type; }
    inline Time LastActionTime() { return LastAction->time; }
    inline Bool LastActionValid() { return LastAction->valid(); }
    inline Bool LastActionValid(Time t) { return LastAction->valid(t); }

    /* have I already called a turn-neck? */
    inline Bool TurnNeckThisCycle() { return TurnNeck.time == CurrentTime ? TRUE : FALSE; }

    Time LastBehaveTime;           // last time into behave
    Time LastActionOpTime;         // last time I could have acted
    Time LastInterruptTime;        // last time into the alarm handler
    int InterruptsThisCycle;       // number of times in so far
    Time LastStartClockTime;       // time when server clock started again
    Time SecondLastStartClockTime; // time when server clock started again

    Bool NewSight;
    Bool NewAction;
    Bool FirstActionOpSinceLastSight;
    Bool SightPredictedEarlyThisCycle();
    Bool GotSightFromCurrentPosition();
    inline Bool TimeToTurnForScan()
    {
        return (SightPredictedEarlyThisCycle() || GotSightFromCurrentPosition()) ? TRUE : FALSE;
    }

    SenseType LastSenseType;

    Vqual ViewQuality;
    Vwidth ViewWidth;
    Vwidth LastViewWidth;
    Time ViewWidthTime;

    Pmode PlayMode;
    Time PlayModeTime;
    Pmode LastPlayMode;
    Kmode KickOffMode;
    int MyScore;
    int TheirScore;

    AngleDeg MyViewAngle(Time time);
    inline AngleDeg MyViewAngle() { return MyViewAngle(CurrentTime); }
    Bool InViewAngle(Time time, AngleDeg ang, float buffer = 5.0);
    inline Bool InViewAngle(AngleDeg ang, float buffer = 5.0) { return InViewAngle(CurrentTime, ang, buffer); }
    int MySightInterval();
    int PredictedNextSightInterval();

    inline void SetMyPos(Vector p, Time t)
    {
        pos = p;
        conf = 1;
        my_pos_time = t;
    }
    inline void SetMyBodyAng(AngleDeg a) { body_ang = a; }
    inline void SetMyNeckRelAng(AngleDeg a) { neck_rel_ang = a; }
    inline void SetMyVel(Vector v, Time t)
    {
        vel = v;
        vel_conf = 1;
        my_vel_time = t;
    }

    void SetMySensedInfo(float st, float e, float sp, float ha, int k, int d, int tu, int sa, int tn, Time ti);
    float GetMySensedSpeed(Time time);
    float GetMySensedStamina(Time time);
    float GetMySensedEffort(Time time);
    float GetMySensedNeckAngle(Time time);
    int GetMySensedKicks(Time time);
    int GetMySensedDashes(Time time);
    int GetMySensedTurns(Time time);
    int GetMySensedSays(Time time);
    int GetMySensedTurnNecks(Time time);

    inline float MyStamina() { return stamina; }
    inline float MyEffort() { return effort; }
    inline float MyRecovery() { return recovery; }
    inline Bool Tired()
    {
        return (MyStamina() < EffortDecThreshold + SP_max_power + CP_tired_buffer) ? TRUE : FALSE;
    }

    float RecoveryDecThreshold;
    float EffortDecThreshold;
    float EffortIncThreshold;
    float CorrectDashPowerForStamina(float dash_power, float stamina, float effort, float recovery);
    inline float CorrectDashPowerForStamina(float dash_power, float stamina)
    {
        return CorrectDashPowerForStamina(dash_power, stamina, MyEffort(), MyRecovery());
    }
    inline float CorrectDashPowerForStamina(float dash_power)
    {
        return CorrectDashPowerForStamina(dash_power, MyStamina(), MyEffort(), MyRecovery());
    }
    void UpdatePredictedStaminaWithDash(float *pStamina, float *pEffort,
                                        float *pRecovery, float dash_power);

    inline AngleDeg EffectiveTurn(AngleDeg actual_turn, float my_speed)
    {
        return 1.0 + actual_turn / (1.0 + SP_inertia_moment * my_speed);
    }
    inline AngleDeg EffectiveTurn(AngleDeg actual_turn)
    {
        return EffectiveTurn(actual_turn, MySpeed());
    }
    inline AngleDeg MaxEffectiveTurn() /* how much we'll actually turn if we try max turn */
    {
        return EffectiveTurn(SP_max_moment);
    }
    inline AngleDeg MaxEffectiveTurn(float my_speed) /* how much we'll actually turn if we try max turn */
    {
        return EffectiveTurn(SP_max_moment, my_speed);
    }

    Bool CanFaceAngleFromNeckWithNeck(AngleDeg ang);
    Bool CanFaceAngleFromBodyWithNeck(AngleDeg ang);
    inline Bool CanFaceGlobalAngleWithNeck(AngleDeg ang)
    {
        return CanFaceAngleFromBodyWithNeck(GetNormalizeAngleDeg(ang - MyBodyAng()));
    }
    inline Bool CanFacePointWithNeck(Vector pt)
    {
        return CanFaceAngleFromBodyWithNeck(AngleToFromBody(pt));
    }

    Bool CanSeeAngleFromNeckWithNeck(AngleDeg ang);
    Bool CanSeeAngleFromBodyWithNeck(AngleDeg ang);
    inline Bool CanSeeGlobalAngleWithNeck(AngleDeg ang)
    {
        return CanSeeAngleFromBodyWithNeck(GetNormalizeAngleDeg(ang - MyBodyAng()));
    }
    inline Bool CanSeePointWithNeck(Vector pt)
    {
        return CanSeeAngleFromBodyWithNeck(AngleToFromBody(pt));
    }

    inline AngleDeg LimitTurnNeckAngle(AngleDeg turn_ang)
    {
        return MinMax(SP_min_neck_angle - MyNeckRelAng(), GetNormalizeAngleDeg(turn_ang),
                      SP_max_neck_angle - MyNeckRelAng());
    }

    inline Bool SensedInfoKnown(Time time) { return (Bool)(sense_time == time || prev_sense_time == time); }

    inline Vector MyPos()
    {
        if (!MyConf())
            my_error("Don't know my pos");
        return pos;
    }
    inline float MyX() { return MyPos().x; }
    inline float MyY() { return MyPos().y; }
    inline AngleDeg MyBodyAng() { return body_ang; }
    inline AngleDeg MyNeckGlobalAng() { return GetNormalizeAngleDeg(MyBodyAng() + MyNeckRelAng()); }
    inline AngleDeg MyNeckRelAng() { return neck_rel_ang; }
    inline Vector MyVel()
    {
        if (!MyVelConf())
            my_error("Don't know my vel");
        return vel;
    }
    inline float MySpeed() { return MyVel().mod(); }
    inline float MyDir() { return MyVel().dir(); }
    inline float MyConf() { return (conf > CP_min_valid_conf) ? conf : 0; }
    inline float MyVelConf() { return (vel_conf > CP_min_valid_conf) ? vel_conf : 0; }
    inline Time MyPosTime() { return my_pos_time; }
    inline Time MyVelTime() { return my_vel_time; }
    inline Time MyUpdateTime() { return Min(my_pos_time, my_vel_time); }

    inline float DistanceTo(Vector vec)
    {
        if (!MyConf())
            my_error("DistanceTo: no conf");
        return (vec - pos).mod();
    }
    inline AngleDeg AngleToFromBody(Vector vec)
    {
        if (!MyConf())
            my_error("AngleTo: no conf");
        AngleDeg angto = (vec - pos).dir() - MyBodyAng();
        NormalizeAngleDeg(&angto);
        return angto;
    }
    inline AngleDeg AngleToFromNeck(Vector vec)
    {
        if (!MyConf())
            my_error("AngleTo: no conf");
        AngleDeg angto = (vec - pos).dir() - MyNeckGlobalAng();
        NormalizeAngleDeg(&angto);
        return angto;
    }
    inline AngleDeg AngleToGlobal(Vector vec)
    {
        if (!MyConf())
            my_error("AngleTo: no conf");
        AngleDeg angto = (vec - pos).dir();
        NormalizeAngleDeg(&angto);
        return angto;
    }
    inline Vector BodyPolar2Gpos(float dist, AngleDeg ang)
    {
        if (!MyConf())
            my_error("Polar2Gpos: no conf");
        Vector rpos = Polar2Vector(dist, ang);
        return MyPos() + rpos.rotate(MyBodyAng());
    }
    inline Vector NeckPolar2Gpos(float dist, AngleDeg ang)
    {
        if (!MyConf())
            my_error("Polar2Gpos: no conf");
        Vector rpos = Polar2Vector(dist, ang);
        return MyPos() + rpos.rotate(MyNeckGlobalAng());
    }

    inline Bool AtPos(Vector p, float buffer)
    {
        return (DistanceTo(p) <= buffer) ? TRUE : FALSE;
    }
    inline Bool AtPos(Vector p)
    {
        return AtPos(p, CP_at_point_buffer);
    }

    inline Vector Global2RelativeToMyBody(Vector v) { return v.Global2Relative(MyPos(), MyBodyAng()); }
    inline Vector Global2RelativeToMyNeck(Vector v) { return v.Global2Relative(MyPos(), MyNeckGlobalAng()); }
    inline Vector RelativeToMyBody2Global(Vector v) { return v.Relative2Global(MyPos(), MyBodyAng()); }
    inline Vector RelativeToMyNeck2Global(Vector v) { return v.Relative2Global(MyPos(), MyNeckGlobalAng()); }

    /* Cone stuff */
    Bool AmIInCone(float wid_dist_ratio, Vector end, Vector vert)
    {
        return IsPointInCone(MyPos(), wid_dist_ratio, end, vert);
    }

    /* Predicted Position */
    Vector MyPredictedPositionWithTurn(float angle,
                                       int steps = 1, float dash_power = 0,
                                       bool with_turn = TRUE,
                                       int idle_cycles = 0);
    inline Vector MyPredictedPosition(int steps = 1, float dash_power = 0,
                                      int idle_cycles = 0)
    {
        return MyPredictedPositionWithTurn(0.0, steps, dash_power, FALSE, idle_cycles);
    }
    Vector MyPredictedPositionWithQueuedActions();
    AngleDeg MyPredictedBodyAngleWithQueuedActions();
    AngleDeg PredictedPointRelAngFromBodyWithQueuedActions(Vector point);

    Vector MyPredictedPositionAtMaxSpeed(int steps = 1);

    /* How long to point */
    int PredictedCyclesToPoint(Vector pt, float dash_power);
    inline int PredictedCyclesToPoint(Vector pt)
    {
        return PredictedCyclesToPoint(pt, SP_max_power);
    }

    /* how long it takes to turn to an angle */
    int NumTurnsToAngle(float targ_body_ang, float curr_body_ang, float curr_speed);
    int NumTurnsToAngle(float targ_body_ang, float curr_body_ang)
    {
        return NumTurnsToAngle(targ_body_ang, curr_body_ang, MySpeed());
    }
    int NumTurnsToAngle(float targ_body_ang)
    {
        return NumTurnsToAngle(targ_body_ang, MyBodyAng(), MySpeed());
    }

    /* Which side of the field am I on? */
    inline Fside LocationSide(Vector p) { return p.y < 0 ? FS_Left : FS_Right; }
    inline Fside MyLocationSide() { return LocationSide(MyPos()); }

    int kicks;
    int dashes;
    int turns;
    int says;
    int turn_necks;

protected:
    Time my_pos_time;
    Time my_vel_time;

    AngleDeg my_last_body_ang;
    AngleDeg my_last_neck_global_ang;
    AngleDeg my_last_neck_rel_ang;

private:
    Vector pos;
    AngleDeg body_ang;
    AngleDeg neck_rel_ang; /* neck angle relative to the body */
    Vector vel;
    float conf;
    float vel_conf;

    float stamina;
    float effort;
    float recovery;

    Time sense_time;
    Time prev_sense_time;

    float last_speed; /* from sense_body */
    float prev_speed;

    float last_stamina;
    float prev_stamina;

    float last_effort;
    float prev_effort;

    float last_neck_rel_ang;
    float prev_neck_rel_ang;

    int last_kicks;
    int prev_kicks;
    int last_dashes;
    int prev_dashes;
    int last_turns;
    int prev_turns;
    int last_says;
    int prev_says;
    int last_turn_necks;
    int prev_turn_necks;

    struct timeval real_time_of_last_send; // the real time of the last send message
};

/* -*- Mode: C++ -*- */
/* MemPosition.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

class Object
{
public:
    /* Methods to get the position  */
    float get_dist(); /* relative */
    AngleDeg get_ang_from_body();
    AngleDeg get_ang_from_neck();

    inline Vector get_abs_pos()
    {
        if (!pos_valid())
            my_error("can't get_abs_pos %d", type);
        return gpos;
    }
    Vector get_rel_to_body_pos();
    Vector get_rel_to_neck_pos();
    inline float get_x() { return get_abs_pos().x; }
    inline float get_y() { return get_abs_pos().y; }
    inline float pos_valid() { return ((gconf >= min_valid_conf) ? gconf : 0); }

    /* Methods to change the position                      */
    /* Changes do NOT take effect until update() is called */
    virtual void set_polar_from_neck(float d, float a, Time time);
    virtual void set_angle_from_neck(AngleDeg a, Time time);
    virtual void set_chinfo(float dist, float dir, Time time);

    /* Method which processes all changes */
    virtual void update();
    virtual void reset(); /* 0 confidences */
    virtual void clear_seen();
    virtual void sanitize_times();

    ObjType type;

    inline Time GetSeenTime() { return seen_time; }
    inline Time GetSeenMovingTime() { return chtime; }

protected:
    float max_conf;
    float min_valid_conf;

    Bool seen;
    Bool seen_moving;
    Time seen_time;

    /* Based on the object's current position, should it be seen? */
    virtual Bool in_view_range(AngleDeg view_ang, float angle_buffer, float distance_buffer);

    float dist; /* polar */
    AngleDeg ang_from_neck;
    Time pdtime;
    Time patime;

    inline Time ptime() { return Min(pdtime, patime); }

    float distch; /* polar */
    float dirch;
    Time chtime;

    Vector rbpos; /* relative to body */
    Time rbtime;
    Vector rnpos; /* relative to neck */
    Time rntime;

    Vector gpos; /* global   */
    Time gtime;
    float gconf;
};

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

class StationaryObject : public Object
{
public:
    void Initialize(MarkerType m, Vector pos, float max, float min_valid, Bool rotate);
    void Initialize(SideLine l, Vector pos, float max, float min_valid, Bool rotate);

    Vector get_my_pos(AngleDeg my_neck_global_ang);
    AngleDeg get_my_neck_global_ang();
    Vector get_my_vel(AngleDeg my_neck_global_ang);

    int object_id;
};

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

class MobileObject : public Object
{
public:
    virtual void Initialize(ObjType t, float max, float min_valid, float decay, float motion, float speed);

    AngleDeg get_rel_to_body_heading();
    AngleDeg get_rel_to_neck_heading();
    Vector get_rel_to_body_vel(); /* relative direction, absolute speed */
    Vector get_rel_to_neck_vel(); /* relative direction, absolute speed */

    /* Methods to get the velocity */
    inline Vector get_abs_vel() /* global   */
    {
        if (!vel_valid())
            my_error("can't get_abs_vel %d", type);
        return gvel;
    }

    inline float get_speed() { return get_abs_vel().mod(); }
    inline AngleDeg get_abs_heading() { return get_abs_vel().dir(); }
    inline float vel_valid() { return ((gvconf >= min_valid_conf) ? gvconf : 0); }
    virtual inline Bool moving() { return vel_valid() && get_speed() != 0 ? TRUE : FALSE; }

    /* Methods to change the position                      */
    /* Changes do NOT take effect until update() is called */
    virtual void set_heard_info(float x, float y, float conf, float dist, Time time);
    virtual void set_heard_info(float x, float y, float pconf, float dx, float dy, float vconf,
                                float dist, Time time);

    /* Method which processes all changes */
    virtual void update_kick(Time time) = 0; /* pure virtual function */
    virtual void estimate_pos(Time time);
    void estimate_vel(Time time);
    virtual void update(Time time);
    virtual void update_seen(Time time);
    void update_estimate(Time time);
    void update_heard(Time time);
    virtual void reset();
    virtual void clear_seen();
    virtual void forget();
    virtual void sanitize_times();

    Vector estimate_future_pos(int steps, Vector extra_vel = 0.0, Vector extra_vel_per = 0.0);

protected:
    float conf_decay;
    float motion_decay;

    float max_speed;

    Vector gvel; /* global   */

    Time gvtime;
    float gvconf;

private:
    Bool heard;
    Bool heard_moving;

    Vector rbvel; /* body -- relative direction, absolute speed */
    Time rbvtime;
    Vector rnvel; /* neck -- relative direction, absolute speed */
    Time rnvtime;

    Vector hpos; /* heard */
    Vector hvel;
    Time hptime;
    Time hvtime;
    float hpdist;
    float hvdist;
    float hpconf;
    float hvconf;
};

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

class BallObject : public MobileObject
{
public:
    void Initialize(float max, float min_valid, float decay, float motion, float max_speed);

    Bool moving();
    Bool kickable(float buffer = 0.0);
    Bool catchable();
    float get_kick_rate(Time time);

    void update_kick(Time time);
    void estimate_pos(Time time);
    void update(Time time);
    void update_seen(Time time);

    float calc_kick_rate(float dist, float ang);
    float calc_kick_rate() { return calc_kick_rate(get_dist(), get_ang_from_body()); }

    void set_past_kick(float pow, AngleDeg ang, Time t);
    void forget_past_kick(Time t);

protected:
    Bool in_view_range(AngleDeg view_ang, float angle_buffer, float distance_buffer);

private:
    float effective_kick_rate;
    Time ektime;

    Bool use_pos_based_vel_estimate;
    Time pos_based_vel_time;

    Time last_seen_time;
    Vector last_seen_pos;

    /* we keep a record of the most recent kicks we try to do so that
       we can estimate the velocity of the ball from subsequent positions */
#define num_past_kicks (7)
    int past_kick_idx;
    struct PastKick
    {
        Time time;
        Vector kick_effect;
    } past_kicks[num_past_kicks];

    int past_kick_inc(int i)
    {
        return (i + 1) % num_past_kicks;
    }
    int past_kick_dec(int i)
    {
        return ((i - 1) + num_past_kicks) % num_past_kicks;
    }
};

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

class PlayerObject : public MobileObject
{
public:
    void Initialize(float max, float min_valid, float decay, float motion, float max_speed);

    char side;
    Unum unum;

    AngleDeg get_rel_to_body_body_ang();
    AngleDeg get_rel_to_body_neck_ang();
    AngleDeg get_rel_to_neck_body_ang();
    AngleDeg get_rel_to_neck_neck_ang();

    AngleDeg get_neck_rel_ang(); /* The player's relative neck angle */

    /* Methods to get the velocity */
    inline AngleDeg get_abs_body_ang() /* global   */
    {
        if (!body_ang_valid())
            my_error("can't get_abs_body_ang %d", type);
        return gbang;
    }
    inline AngleDeg get_abs_neck_ang() /* global   */
    {
        if (!neck_ang_valid())
            my_error("can't get_abs_neck_ang %d", type);
        return gnang;
    }
    inline float body_ang_valid() { return ((gbaconf >= min_valid_conf) ? gbaconf : 0); }
    inline float neck_ang_valid() { return ((gnaconf >= min_valid_conf) ? gnaconf : 0); }

    /* Methods to change the position                      */
    /* Changes do NOT take effect until update() is called */
    void set_body_ang_from_neck(AngleDeg body_ang, Time time);
    void set_neck_ang_from_neck(AngleDeg neck_ang, Time time);
    inline void set_body_and_neck_ang_from_neck(AngleDeg body_ang, AngleDeg neck_ang, Time time)
    {
        set_body_ang_from_neck(body_ang, time);
        set_neck_ang_from_neck(neck_ang, time);
    }

    void set_heard_info_w_angs(float x, float y, float pconf, float dx, float dy, float vconf,
                               AngleDeg bang, float bconf, AngleDeg nang, float nconf,
                               float dist, Time time);

    /* Method which processes all changes */
    void update(Time time);
    void update_seen_body_and_neck_angs(Time time);
    void update_estimate(Time time);
    void update_heard(Time time);
    void reset();
    void clear_seen();
    void forget();
    void sanitize_times();

    /* This is a pure virtual MobileObject fn.  Shouldn't be called for playerobject */
    inline void update_kick(Time time) { my_error("Can't kick players %d", time.t); }

protected:
    Bool in_view_range(AngleDeg view_ang, float angle_buffer, float distance_buffer);

private:
    Bool seen_body_ang;
    Bool heard_body_ang;
    Bool seen_neck_ang;
    Bool heard_neck_ang;

    AngleDeg rnbang; /* relative to neck body angle */
    Time rnbatime;
    AngleDeg rnnang; /* relative to neck neck angle */
    Time rnnatime;

    AngleDeg gbang; /* global body angle */
    float gbaconf;
    Time gbatime;
    AngleDeg gnang; /* global neck angle */
    float gnaconf;
    Time gnatime;

    AngleDeg hbang; /* heard body angle (global) */
    float hbaconf;
    float hbadist;
    Time hbatime;
    AngleDeg hnang; /* heard neck angle (global) */
    float hnaconf;
    float hnadist;
    Time hnatime;
};

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

class TempPlayerObject
{
public:
    float dist; /* relative */
    AngleDeg ang_from_neck;

    Time time;
    char side;

    inline void set(char s, float d, float a, Time t)
    {
        time = t;
        side = s;
        dist = d;
        ang_from_neck = a;
    }
};

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

#define MAX_MARKERS 100
#define MAX_PLAYERS 25

#define Unum_Unknown 0
#define Unum_Teamless 100

class PositionInfo : public PlayerInfo
{
public:
    ~PositionInfo();
    void Initialize();

    void SeeLine(SideLine l, float dist, float ang, Time tm);
    void SeeLine(SideLine l, float ang, Time tm);

    void SeeMarker(MarkerType marker, float dist, float ang, Time tm);
    void SeeMarker(MarkerType marker, float ang, Time tm);
    void SeeMarker(MarkerType marker, float dist, float ang, float distChng, float dirChng, Time tm);

    void SeeBall(float ang, Time tm);
    void SeeBall(float dist, float ang, Time tm);
    void SeeBall(float dist, float ang, float distChng, float dirChng, Time tm);

    void SeePlayer(float ang, Time time);
    void SeePlayer(float dist, float ang, Time time);
    void SeePlayer(char side, float ang, Time time);
    void SeePlayer(char side, float dist, float ang, Time time);
    void SeePlayer(char side, Unum num, float dist, float ang,
                   float distChng, float dirChng, float bodydir, float neckdir, Time time);
    void SeePlayer(char side, Unum num, float dist, float ang, Time time);
    void SeePlayer(char side, Unum num, float ang, Time time);

    void HearBall(float x, float y, float conf, float dist, Time time);
    void HearBall(float x, float y, float pconf, float dx, float dy, float vconf, float dist, Time time);

    void HearPlayer(char side, Unum num, float x, float y, float conf, float dist, Time time);
    void HearPlayer(char side, Unum num,
                    float x, float y, float pconf, float dx, float dy, float vconf, float dist, Time time);
    void HearPlayer(char side, Unum num, float x, float y, float pconf, float dx, float dy, float vconf,
                    AngleDeg bang, float bconf, AngleDeg nang, float nconf, float dist, Time time);

    /* Access shortcuts -- Markers */
    inline float MarkerDistance(MarkerType m) { return GetMarker(m)->get_dist(); }
    inline float MarkerDistance2(MarkerType m) { return Sqr(GetMarker(m)->get_dist()); }
    inline AngleDeg MarkerAngleFromBody(MarkerType m) { return GetMarker(m)->get_ang_from_body(); }
    inline AngleDeg MarkerAngleFromNeck(MarkerType m) { return GetMarker(m)->get_ang_from_neck(); }
    inline Vector MarkerPosition(MarkerType m) { return GetMarker(m)->get_abs_pos(); }
    inline float MarkerX(MarkerType m) { return GetMarker(m)->get_x(); }
    inline float MarkerY(MarkerType m) { return GetMarker(m)->get_y(); }
    inline float MarkerPositionValid(MarkerType m) { return GetMarker(m)->pos_valid(); }

    /* Access shortcuts -- Ball */
    inline float BallDistance() { return GetBall()->get_dist(); }
    inline AngleDeg BallAngleFromBody() { return GetBall()->get_ang_from_body(); }
    inline AngleDeg BallAngleFromNeck() { return GetBall()->get_ang_from_neck(); }
    inline Vector BallAbsolutePosition() { return GetBall()->get_abs_pos(); }
    inline Vector BallRelativeToBodyPosition() { return GetBall()->get_rel_to_body_pos(); }
    inline Vector BallRelativeToNeckPosition() { return GetBall()->get_rel_to_neck_pos(); }
    inline float BallX() { return GetBall()->get_x(); }
    inline float BallY() { return GetBall()->get_y(); }
    inline float BallPositionValid() { return GetBall()->pos_valid(); }

    inline float BallSpeed() { return GetBall()->get_speed(); }
    inline AngleDeg BallRelativeToBodyHeading() { return GetBall()->get_rel_to_body_heading(); }
    inline AngleDeg BallRelativeToNeckHeading() { return GetBall()->get_rel_to_body_heading(); }
    inline AngleDeg BallAbsoluteHeading() { return GetBall()->get_abs_heading(); }
    inline Vector BallRelativeToBodyVelocity() { return GetBall()->get_rel_to_body_vel(); }
    inline Vector BallRelativeToNeckVelocity() { return GetBall()->get_rel_to_neck_vel(); }
    inline Vector BallAbsoluteVelocity() { return GetBall()->get_abs_vel(); }
    inline float BallVelocityValid() { return GetBall()->vel_valid(); }
    inline Bool BallMoving() { return GetBall()->moving(); }
    inline Bool BallKickable(float buffer = 0.0)
    {
        return GetBall()->kickable(buffer);
    }
    inline Bool BallCatchable()
    {
        return (GetBall()->catchable() && InOwnPenaltyArea() &&
                CP_goalie && PlayMode == PM_Play_On)
                   ? TRUE
                   : FALSE;
    }
    inline float BallKickRate() { return GetBall()->get_kick_rate(CurrentTime); }

    /* Access shortcuts -- Players */
    inline float PlayerDistance(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_dist();
    }
    inline AngleDeg PlayerAngleFromBody(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_ang_from_body();
    }
    inline AngleDeg PlayerAngleFromNeck(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_ang_from_neck();
    }
    inline Vector PlayerAbsolutePosition(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? MyPos() : GetPlayer(s, n)->get_abs_pos();
    }
    inline Vector PlayerRelativeToBodyPosition(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? Vector(0, 0) : GetPlayer(s, n)->get_rel_to_body_pos();
    }
    inline Vector PlayerRelativeToNeckPosition(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? Vector(0, 0) : GetPlayer(s, n)->get_rel_to_neck_pos();
    }
    inline float PlayerX(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? MyX() : GetPlayer(s, n)->get_x();
    }
    inline float PlayerY(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? MyY() : GetPlayer(s, n)->get_y();
    }
    float PlayerPositionValid(char s, Unum n);

    inline float PlayerSpeed(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? MySpeed() : GetPlayer(s, n)->get_speed();
    }
    inline AngleDeg PlayerRelativeToBodyHeading(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_rel_to_body_heading();
    }
    inline AngleDeg PlayerRelativeToNeckHeading(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_rel_to_neck_heading();
    }
    inline AngleDeg PlayerAbsoluteHeading(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? MyDir() : GetPlayer(s, n)->get_abs_heading();
    }
    inline Vector PlayerAbsoluteVelocity(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? MyVel() : GetPlayer(s, n)->get_abs_vel();
    }
    inline Vector PlayerRelativeToBodyVelocity(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? Polar2Vector(MySpeed(), 0) : GetPlayer(s, n)->get_rel_to_body_vel();
    }
    inline Vector PlayerRelativeToNeckVelocity(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? Polar2Vector(MySpeed(), 0) : GetPlayer(s, n)->get_rel_to_neck_vel();
    }
    float PlayerVelocityValid(char s, Unum n);
    inline AngleDeg PlayerRelativeToBodyBodyAngle(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_rel_to_body_body_ang();
    }
    inline AngleDeg PlayerRelativeToNeckBodyAngle(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_rel_to_neck_body_ang();
    }
    inline AngleDeg PlayerRelativeToBodyNeckAngle(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_rel_to_body_neck_ang();
    }
    inline AngleDeg PlayerRelativeToNeckNeckAngle(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_rel_to_neck_neck_ang();
    }
    inline AngleDeg PlayerRelativeNeckAngle(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? 0 : GetPlayer(s, n)->get_neck_rel_ang();
    }
    inline AngleDeg PlayerAbsoluteBodyAngle(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? MyBodyAng() : GetPlayer(s, n)->get_abs_body_ang();
    }
    float PlayerBodyAngleValid(char s, Unum n);
    inline AngleDeg PlayerAbsoluteNeckAngle(char s, Unum n)
    {
        return (s == MySide && n == MyNumber) ? MyNeckGlobalAng() : GetPlayer(s, n)->get_abs_neck_ang();
    }
    float PlayerNeckAngleValid(char s, Unum n);

    /* Access shortcuts -- Teammates */
    inline float TeammateDistance(Unum n) { return PlayerDistance(MySide, n); }
    inline AngleDeg TeammateAngleFromBody(Unum n) { return PlayerAngleFromBody(MySide, n); }
    inline AngleDeg TeammateAngleFromNeck(Unum n) { return PlayerAngleFromNeck(MySide, n); }
    inline Vector TeammateAbsolutePosition(Unum n) { return PlayerAbsolutePosition(MySide, n); }
    inline Vector TeammateRelativeToBodyPosition(Unum n) { return PlayerRelativeToBodyPosition(MySide, n); }
    inline Vector TeammateRelativeToNeckPosition(Unum n) { return PlayerRelativeToNeckPosition(MySide, n); }
    inline float TeammateX(Unum n) { return PlayerX(MySide, n); }
    inline float TeammateY(Unum n) { return PlayerY(MySide, n); }
    inline float TeammatePositionValid(Unum n) { return PlayerPositionValid(MySide, n); }

    inline float TeammateSpeed(Unum n) { return PlayerSpeed(MySide, n); }
    inline AngleDeg TeammateRelativeToBodyHeading(Unum n) { return PlayerRelativeToBodyHeading(MySide, n); }
    inline AngleDeg TeammateRelativeToNeckHeading(Unum n) { return PlayerRelativeToNeckHeading(MySide, n); }
    inline AngleDeg TeammateAbsoluteHeading(Unum n) { return PlayerAbsoluteHeading(MySide, n); }
    inline Vector TeammateAbsoluteVelocity(Unum n) { return PlayerAbsoluteVelocity(MySide, n); }
    inline Vector TeammateRelativeToBodyVelocity(Unum n) { return PlayerRelativeToBodyVelocity(MySide, n); }
    inline Vector TeammateRelativeToNeckVelocity(Unum n) { return PlayerRelativeToNeckVelocity(MySide, n); }
    inline float TeammateVelocityValid(Unum n) { return PlayerVelocityValid(MySide, n); }
    inline AngleDeg TeammateRelativeToBodyBodyAngle(Unum n) { return PlayerRelativeToBodyBodyAngle(MySide, n); }
    inline AngleDeg TeammateRelativeToNeckBodyAngle(Unum n) { return PlayerRelativeToNeckBodyAngle(MySide, n); }
    inline AngleDeg TeammateRelativeToBodyNeckAngle(Unum n) { return PlayerRelativeToBodyNeckAngle(MySide, n); }
    inline AngleDeg TeammateRelativeToNeckNeckAngle(Unum n) { return PlayerRelativeToNeckNeckAngle(MySide, n); }
    inline AngleDeg TeammateRelativeNeckAngle(Unum n) { return PlayerRelativeNeckAngle(MySide, n); }
    inline AngleDeg TeammateAbsoluteBodyAngle(Unum n) { return PlayerAbsoluteBodyAngle(MySide, n); }
    inline float TeammateBodyAngleValid(Unum n) { return PlayerBodyAngleValid(MySide, n); }
    inline AngleDeg TeammateAbsoluteNeckAngle(Unum n) { return PlayerAbsoluteNeckAngle(MySide, n); }
    inline float TeammateNeckAngleValid(Unum n) { return PlayerNeckAngleValid(MySide, n); }

    /* Access shortcuts -- Opponents */
    inline float OpponentDistance(Unum n) { return PlayerDistance(TheirSide, n); }
    inline AngleDeg OpponentAngleFromBody(Unum n) { return PlayerAngleFromBody(TheirSide, n); }
    inline AngleDeg OpponentAngleFromNeck(Unum n) { return PlayerAngleFromNeck(TheirSide, n); }
    inline Vector OpponentAbsolutePosition(Unum n) { return PlayerAbsolutePosition(TheirSide, n); }
    inline Vector OpponentRelativeToBodyPosition(Unum n) { return PlayerRelativeToBodyPosition(TheirSide, n); }
    inline Vector OpponentRelativeToNeckPosition(Unum n) { return PlayerRelativeToNeckPosition(TheirSide, n); }
    inline float OpponentX(Unum n) { return PlayerX(TheirSide, n); }
    inline float OpponentY(Unum n) { return PlayerY(TheirSide, n); }
    inline float OpponentPositionValid(Unum n) { return PlayerPositionValid(TheirSide, n); }

    inline float OpponentSpeed(Unum n) { return PlayerSpeed(TheirSide, n); }
    inline AngleDeg OpponentRelativeToBodyHeading(Unum n) { return PlayerRelativeToBodyHeading(TheirSide, n); }
    inline AngleDeg OpponentRelativeToNeckHeading(Unum n) { return PlayerRelativeToNeckHeading(TheirSide, n); }
    inline AngleDeg OpponentAbsoluteHeading(Unum n) { return PlayerAbsoluteHeading(TheirSide, n); }
    inline Vector OpponentAbsoluteVelocity(Unum n) { return PlayerAbsoluteVelocity(TheirSide, n); }
    inline Vector OpponentRelativeToBodyVelocity(Unum n) { return PlayerRelativeToBodyVelocity(TheirSide, n); }
    inline Vector OpponentRelativeToNeckVelocity(Unum n) { return PlayerRelativeToNeckVelocity(TheirSide, n); }
    inline float OpponentVelocityValid(Unum n) { return PlayerVelocityValid(TheirSide, n); }
    inline AngleDeg OpponentRelativeToBodyBodyAngle(Unum n) { return PlayerRelativeToBodyBodyAngle(TheirSide, n); }
    inline AngleDeg OpponentRelativeToNeckBodyAngle(Unum n) { return PlayerRelativeToNeckBodyAngle(TheirSide, n); }
    inline AngleDeg OpponentRelativeToBodyNeckAngle(Unum n) { return PlayerRelativeToBodyNeckAngle(TheirSide, n); }
    inline AngleDeg OpponentRelativeToNeckNeckAngle(Unum n) { return PlayerRelativeToNeckNeckAngle(TheirSide, n); }
    inline AngleDeg OpponentRelativeNeckAngle(Unum n) { return PlayerRelativeNeckAngle(TheirSide, n); }
    inline AngleDeg OpponentAbsoluteBodyAngle(Unum n) { return PlayerAbsoluteBodyAngle(TheirSide, n); }
    inline float OpponentBodyAngleValid(Unum n) { return PlayerBodyAngleValid(TheirSide, n); }
    inline AngleDeg OpponentAbsoluteNeckAngle(Unum n) { return PlayerAbsoluteNeckAngle(TheirSide, n); }
    inline float OpponentNeckAngleValid(Unum n) { return PlayerNeckAngleValid(TheirSide, n); }

    /* more complex shortcuts */
    /* kickable for other players */
    Bool BallKickableForPlayer(char s, Unum n, float buffer = 0);
    inline Bool BallKickableForTeammate(Unum n, float buffer = 0)
    {
        return BallKickableForPlayer(MySide, n, buffer);
    }
    inline Bool BallKickableForOpponent(Unum n, float buffer = 0)
    {
        return BallKickableForPlayer(TheirSide, n, buffer);
    }

    inline Bool CanFaceBallWithNeck()
    {
        return CanFaceAngleFromBodyWithNeck(BallAngleFromBody());
    }
    inline Bool CanFacePlayerWithNeck(char s, Unum n)
    {
        return CanFaceAngleFromBodyWithNeck(PlayerAngleFromBody(s, n));
    }
    inline Bool CanFaceTeammateWithNeck(Unum n)
    {
        return CanFaceAngleFromBodyWithNeck(TeammateAngleFromBody(n));
    }
    inline Bool CanFaceOpponentWithNeck(Unum n)
    {
        return CanFaceAngleFromBodyWithNeck(OpponentAngleFromBody(n));
    }

    inline Bool CanSeeBallWithNeck()
    {
        return CanSeeAngleFromBodyWithNeck(BallAngleFromBody());
    }
    inline Bool CanSeePlayerWithNeck(char s, Unum n)
    {
        return CanSeeAngleFromBodyWithNeck(PlayerAngleFromBody(s, n));
    }
    inline Bool CanSeeTeammateWithNeck(Unum n)
    {
        return CanSeeAngleFromBodyWithNeck(TeammateAngleFromBody(n));
    }
    inline Bool CanSeeOpponentWithNeck(Unum n)
    {
        return CanSeeAngleFromBodyWithNeck(OpponentAngleFromBody(n));
    }

    /* is player behind us */
    /* three forms for third arg:
       angle: absolute angle indicating the front
       Vector: position whose direction indicates the front
       none: use MyAng */
    /* is point behind us */
    /* three forms for third arg:
       angle: absolute angle indicating the front
       Vector: position whose direction indicates the front
       none: use MyAng */
    inline Bool IsPointBehind(Vector pt, AngleDeg ang)
    {
        return (fabs(GetNormalizeAngleDeg(AngleToFromBody(pt) + MyBodyAng() - ang)) > CP_behind_angle)
                   ? TRUE
                   : FALSE;
    }
    inline Bool IsPointBehind(Vector pt, Vector targPos)
    {
        return IsPointBehind(pt, (targPos - MyPos()).dir());
    }
    inline Bool IsPointBehindBody(Vector pt)
    {
        return (fabs(AngleToFromBody(pt)) > CP_behind_angle)
                   ? TRUE
                   : FALSE;
    }
    inline Bool IsPointBehindNeck(Vector pt)
    {
        return (fabs(AngleToFromNeck(pt)) > CP_behind_angle)
                   ? TRUE
                   : FALSE;
    }

    inline Bool IsPlayerBehind(char s, Unum n, AngleDeg ang)
    {
        if (!PlayerPositionValid(s, n))
            my_error("Can't tell if player is behind if I don't know where he is");
        return (fabs(PlayerAngleFromBody(s, n) + MyBodyAng() - ang) > CP_behind_angle)
                   ? TRUE
                   : FALSE;
    }
    inline Bool IsPlayerBehind(char s, Unum n, Vector targPos)
    {
        return IsPlayerBehind(s, n, (targPos - MyPos()).dir());
    }
    inline Bool IsPlayerBehindFromBody(char s, Unum n)
    {
        if (!PlayerPositionValid(s, n))
            my_error("Can't tell if player is behind if I don't know where he is");
        return (fabs(PlayerAngleFromBody(s, n)) > CP_behind_angle)
                   ? TRUE
                   : FALSE;
    }

    inline Bool IsPlayerBehindFromNeck(char s, Unum n)
    {
        if (!PlayerPositionValid(s, n))
            my_error("Can't tell if player is behind if I don't know where he is");
        return (fabs(PlayerAngleFromNeck(s, n)) > CP_behind_angle)
                   ? TRUE
                   : FALSE;
    }

    inline Bool IsTeammateBehind(Unum n, AngleDeg ang)
    {
        return IsPlayerBehind(MySide, n, ang);
    }
    inline Bool IsTeammateBehind(Unum n, Vector targPos)
    {
        return IsPlayerBehind(MySide, n, targPos);
    }
    inline Bool IsTeammateBehindFromBody(Unum n)
    {
        return IsPlayerBehindFromBody(MySide, n);
    }
    inline Bool IsTeammateBehindFromNeck(Unum n)
    {
        return IsPlayerBehindFromNeck(MySide, n);
    }

    inline Bool IsOpponentBehind(Unum n, AngleDeg ang)
    {
        return IsPlayerBehind(TheirSide, n, ang);
    }
    inline Bool IsOpponentBehind(Unum n, Vector targPos)
    {
        return IsPlayerBehind(TheirSide, n, targPos);
    }
    inline Bool IsOpponentBehindFromBody(Unum n)
    {
        return IsPlayerBehindFromBody(TheirSide, n);
    }
    inline Bool IsOpponentBehindFromNeck(Unum n)
    {
        return IsPlayerBehindFromNeck(TheirSide, n);
    }

    /* Relative markers */
    MarkerType RM_My_Goal;
    MarkerType RM_Their_Goal;
    MarkerType RM_LB_Flag;
    MarkerType RM_LC_Flag;
    MarkerType RM_LF_Flag;
    MarkerType RM_RB_Flag;
    MarkerType RM_RC_Flag;
    MarkerType RM_RF_Flag;
    MarkerType RM_My_PC_Flag;    /* Center of my penalty area */
    MarkerType RM_Their_PC_Flag; /* Center of theirs          */

    /* predicted positions */
    Vector BallPredictedPosition(int steps, float kick_power, AngleDeg kick_dir);
    Vector BallPredictedPositionWithQueuedActions(int steps = 1);
    inline Vector BallPredictedPosition(int steps = 1)
    {
        return GetBall()->estimate_future_pos(steps);
    }

    inline Vector PlayerPredictedPosition(char side, Unum num, int steps = 1, Vector dash_per = Vector(0, 0))
    {
        return GetPlayer(side, num)->estimate_future_pos(steps, 0, dash_per);
    }
    inline Vector TeammatePredictedPosition(Unum num, int steps = 1, Vector dash_per = Vector(0, 0))
    {
        return PlayerPredictedPosition(MySide, num, steps, dash_per);
    }
    inline Vector OpponentPredictedPosition(Unum num, int steps = 1, Vector dash_per = Vector(0, 0))
    {
        return PlayerPredictedPosition(TheirSide, num, steps, dash_per);
    }

    Bool BallWillBeKickable(int steps = 1, float dash_power = 0, float buffer = 0);

    Bool WillKickBeCollision(float kick_power, AngleDeg kick_dir, float buffer = 0)
    {
        return ((BallPredictedPosition(1, kick_power, kick_dir) - MyPredictedPosition()).mod() <=
                SP_player_size + SP_ball_size + buffer)
                   ? TRUE
                   : FALSE;
    }
    Bool WillDashBeCollision(float dash_power, float buffer = 0)
    {
        return ((BallPredictedPosition() - MyPredictedPosition(1, dash_power)).mod() <=
                SP_player_size + SP_ball_size + buffer)
                   ? TRUE
                   : FALSE;
    }

    /* PredictedCyclesToPoint */
    /* a buffer of 0 is a bad idea! */
    int PlayerPredictedCyclesToPoint(char side, Unum num, Vector pt,
                                     float dash_power, float buffer);
    int PlayerPredictedCyclesToPoint(char side, Unum num, Vector pt, float dash_power)
    {
        return PlayerPredictedCyclesToPoint(side, num, pt, dash_power, CP_at_point_buffer);
    }
    int PlayerPredictedCyclesToPoint(char side, Unum num, Vector pt)
    {
        return PlayerPredictedCyclesToPoint(side, num, pt, SP_max_power);
    }
    int TeammatePredictedCyclesToPoint(Unum num, Vector pt, float dash_power, float buffer)
    {
        return PlayerPredictedCyclesToPoint(MySide, num, pt, dash_power, buffer);
    }
    int TeammatePredictedCyclesToPoint(Unum num, Vector pt, float dash_power)
    {
        return PlayerPredictedCyclesToPoint(MySide, num, pt, dash_power);
    }
    int TeammatePredictedCyclesToPoint(Unum num, Vector pt)
    {
        return PlayerPredictedCyclesToPoint(MySide, num, pt);
    }
    int OpponentPredictedCyclesToPoint(Unum num, Vector pt, float dash_power, float buffer)
    {
        return PlayerPredictedCyclesToPoint(TheirSide, num, pt, dash_power, buffer);
    }
    int OpponentPredictedCyclesToPoint(Unum num, Vector pt, float dash_power)
    {
        return PlayerPredictedCyclesToPoint(TheirSide, num, pt, dash_power);
    }
    int OpponentPredictedCyclesToPoint(Unum num, Vector pt)
    {
        return PlayerPredictedCyclesToPoint(TheirSide, num, pt);
    }

    Unum PlayerWithBall(float buffer = 0); /* with ball means actually kickable      */
    Unum TeammateWithBall(float buffer = 0);
    Unum OpponentWithBall(float buffer = 0);
    char TeamWithBall(float buffer = 0);

    /* Direct access functions */
    inline StationaryObject *GetMarker(MarkerType m) { return &Marker[m]; }

    inline BallObject *GetBall() { return &Ball; }

    inline PlayerObject **Teammates() { return MyPlayer; }
    inline PlayerObject **Opponents() { return TheirPlayer; }
    inline PlayerObject **TeamlessPlayers() { return TeamlessPlayer; }

    inline int NumTeammates() { return num_my_players; }
    inline int NumOpponents() { return num_their_players; }
    inline int NumTeamlessPlayers() { return num_teamless_players; }
    inline int NumPlayers() { return num_my_players + num_their_players + num_teamless_players; }

    PlayerObject *GetTeammate(Unum num);
    PlayerObject *GetOpponent(Unum num);
    PlayerObject *GetPlayer(char side, Unum num);

    /* Player data-structure management */
    PlayerObject *GetNewPlayer(char side, Unum num);
    Bool ForgetAPlayer(char side);
    void CleanMyPlayers();
    void CleanTheirPlayers();
    void CleanTeamlessPlayers();
    void CleanAllPlayers();
    Bool ResetMyDuplicatePlayers();
    Bool ResetTheirDuplicatePlayers();
    void ClearSeenInfo();

    MarkerType ClosestMarker;
    MarkerType ClosestMotionMarker;

    MarkerType ClosestGoal();
    MarkerType ClosestFlagTo();

    int SortPlayersBy(char side, char keyFunc, float keyNum, Unum *players);
    /* side = m,t,b (both) keyFunc = d,a
(dist/ang), keyNum = central val.
players = array.  Return num found */
    int SortPlayersByDistanceToPoint(char side, Vector point, Unum *players);
    int SortPlayersByDistanceToLine(char side, Line line, Unum *players,
                                    Bool TestEndPoints = FALSE, Vector ep1 = 0, Vector ep2 = 0);

    int NumTeammatesWithin(float Dist, Vector of_pos);
    inline int NumTeammatesWithin(float Dist)
    {
        if (!MyConf())
            my_error("nope");
        return NumTeammatesWithin(Dist, MyPos());
    }
    int NumOpponentsWithin(float Dist, Vector of_pos);
    inline int NumOpponentsWithin(float Dist)
    {
        if (!MyConf())
            my_error("nope");
        return NumOpponentsWithin(Dist, MyPos());
    }
    inline int NumPlayersWithin(float Dist, Vector of_pos)
    {
        return NumTeammatesWithin(Dist, of_pos) + NumOpponentsWithin(Dist, of_pos);
    }
    inline int NumPlayersWithin(float Dist)
    {
        if (!MyConf())
            my_error("nope");
        return NumPlayersWithin(Dist, MyPos());
    }

    /* These are needed for the Decision Tree, Dodge player -- use dist and ang */
    int NumTeammatesWithin(float Dist, AngleDeg Ang, float ofDist, AngleDeg ofAng);
    inline int NumTeammatesWithin(float Dist, AngleDeg Ang) { return NumTeammatesWithin(Dist, Ang, 0, 0); }
    int NumOpponentsWithin(float Dist, AngleDeg Ang, float ofDist, AngleDeg ofAng);
    inline int NumOpponentsWithin(float Dist, AngleDeg Ang) { return NumOpponentsWithin(Dist, Ang, 0, 0); }
    inline int NumPlayersWithin(float Dist, AngleDeg Ang, float ofDist, AngleDeg ofAng)
    {
        return NumTeammatesWithin(Dist, Ang, ofDist, ofAng) + NumOpponentsWithin(Dist, Ang, ofDist, ofAng);
    }
    inline int NumPlayersWithin(float Dist, AngleDeg Ang) { return NumPlayersWithin(Dist, Ang, 0, 0); }

    PlayerObject *GetPlayerWithin(float Dist, Vector ofPos);
    PlayerObject *GetPlayerWithin(float Dist, AngleDeg Ang, float ofDist, AngleDeg ofAng);

    int NumOpponentsInCone(float wid_dist_ratio, Vector end, Vector vert);
    int NumTeammatesInCone(float wid_dist_ratio, Vector end, Vector vert, Bool IncludeMe = FALSE);
    inline int NumPlayersInCone(float wid_dist_ratio, Vector end, Vector vert)
    {
        return NumOpponentsInCone(wid_dist_ratio, end, vert) +
               NumTeammatesInCone(wid_dist_ratio, end, vert);
    }

    int NumOpponentsInCone(float wid_dist_ratio, Vector end)
    {
        return NumOpponentsInCone(wid_dist_ratio, end, MyPos());
    }
    int NumTeammatesInCone(float wid_dist_ratio, Vector end)
    {
        return NumTeammatesInCone(wid_dist_ratio, end, MyPos());
    }
    int NumPlayersInCone(float wid_dist_ratio, Vector end)
    {
        return NumPlayersInCone(wid_dist_ratio, end, MyPos());
    }

    /* which is 'm' (my team), 't' (thier team), 'b' (both) */
    /* this should probably be private */
    int NumPlayersInConeToPlayer(char which,
                                 float wid_dist_ratio, char side, Unum num,
                                 float extra_len, Vector vert);

    int NumPlayersInConeToPlayer(float wid_dist_ratio, char side, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('b', wid_dist_ratio, side, num, extra_len, vert);
    }
    int NumTeammatesInConeToPlayer(float wid_dist_ratio, char side, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('m', wid_dist_ratio, side, num, extra_len, vert);
    }
    int NumOpponentsInConeToPlayer(float wid_dist_ratio, char side, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('t', wid_dist_ratio, side, num, extra_len, vert);
    }
    int NumPlayersInConeToTeammate(float wid_dist_ratio, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('b', wid_dist_ratio, MySide, num, extra_len, vert);
    }
    int NumTeammatesInConeToTeammate(float wid_dist_ratio, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('m', wid_dist_ratio, MySide, num, extra_len, vert);
    }
    int NumOpponentsInConeToTeammate(float wid_dist_ratio, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('t', wid_dist_ratio, MySide, num, extra_len, vert);
    }
    int NumPlayersInConeToOpponent(float wid_dist_ratio, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('b', wid_dist_ratio, TheirSide, num, extra_len, vert);
    }
    int NumTeammatesInConeToOpponent(float wid_dist_ratio, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('m', wid_dist_ratio, TheirSide, num, extra_len, vert);
    }
    int NumOpponentsInConeToOpponent(float wid_dist_ratio, Unum num, float extra_len, Vector vert)
    {
        return NumPlayersInConeToPlayer('t', wid_dist_ratio, TheirSide, num, extra_len, vert);
    }

    int NumPlayersInConeToPlayer(float wid_dist_ratio, char side, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('b', wid_dist_ratio, side, num, extra_len, MyPos());
    }
    int NumTeammatesInConeToPlayer(float wid_dist_ratio, char side, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('m', wid_dist_ratio, side, num, extra_len, MyPos());
    }
    int NumOpponentsInConeToPlayer(float wid_dist_ratio, char side, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('t', wid_dist_ratio, side, num, extra_len, MyPos());
    }
    int NumPlayersInConeToTeammate(float wid_dist_ratio, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('b', wid_dist_ratio, MySide, num, extra_len, MyPos());
    }
    int NumTeammatesInConeToTeammate(float wid_dist_ratio, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('m', wid_dist_ratio, MySide, num, extra_len, MyPos());
    }
    int NumOpponentsInConeToTeammate(float wid_dist_ratio, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('t', wid_dist_ratio, MySide, num, extra_len, MyPos());
    }
    int NumPlayersInConeToOpponent(float wid_dist_ratio, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('b', wid_dist_ratio, TheirSide, num, extra_len, MyPos());
    }
    int NumTeammatesInConeToOpponent(float wid_dist_ratio, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('m', wid_dist_ratio, TheirSide, num, extra_len, MyPos());
    }
    int NumOpponentsInConeToOpponent(float wid_dist_ratio, Unum num, float extra_len = 0.0)
    {
        return NumPlayersInConeToPlayer('t', wid_dist_ratio, TheirSide, num, extra_len, MyPos());
    }

    Unum ClosestTeammateTo(Vector p, Bool include_me = TRUE);
    Unum ClosestOpponentTo(Vector p);
    Vector ClosestTeamlessPlayerPosition();
    Unum ClosestTeammateToBall(Bool include_me = TRUE);
    Unum ClosestOpponentToBall();
    float ClosestTeammateToBallDistance(Bool include_me = TRUE);

    inline int NumTeammatesCloserTo(Vector p) { return NumTeammatesWithin(DistanceTo(p), p) - 1; }
    inline int NumOpponentsCloserTo(Vector p) { return NumOpponentsWithin(DistanceTo(p), p); }
    inline int NumPlayersCloserTo(Vector p) { return NumTeammatesCloserTo(p) + NumOpponentsCloserTo(p); }
    inline int NumTeammatesCloserToBall() { return NumTeammatesCloserTo(BallAbsolutePosition()); }
    inline int NumOpponentsCloserToBall() { return NumOpponentsCloserTo(BallAbsolutePosition()); }
    inline int NumPlayersCloserToBall() { return NumPlayersCloserTo(BallAbsolutePosition()); }

    float PlayerDistanceTo(char s, Unum n, Vector p);
    inline float TeammateDistanceTo(Unum n, Vector p) { return PlayerDistanceTo(MySide, n, p); }
    inline float OpponentDistanceTo(Unum n, Vector p) { return PlayerDistanceTo(TheirSide, n, p); }

    float PlayerDistance2To(char s, Unum n, Vector p);
    inline float TeammateDistance2To(Unum n, Vector p) { return PlayerDistance2To(MySide, n, p); }
    inline float OpponentDistance2To(Unum n, Vector p) { return PlayerDistance2To(TheirSide, n, p); }

    float PlayerDistanceToLine(char s, Unum n, Line l);
    inline float TeammateDistanceToLine(Unum n, Line l) { return PlayerDistanceToLine(MySide, n, l); }
    inline float OpponentDistanceToLine(Unum n, Line l) { return PlayerDistanceToLine(TheirSide, n, l); }

    float PlayerDistance2ToLine(char s, Unum n, Line l);
    inline float TeammateDistance2ToLine(Unum n, Line l) { return PlayerDistance2ToLine(MySide, n, l); }
    inline float OpponentDistance2ToLine(Unum n, Line l) { return PlayerDistance2ToLine(TheirSide, n, l); }

    inline float PlayerDistanceToBall(char side, Unum num)
    {
        if (!BallPositionValid())
            my_error("Ball unknown");
        return PlayerDistanceTo(side, num, BallAbsolutePosition());
    }
    inline float TeammateDistanceToBall(Unum num) { return PlayerDistanceToBall(MySide, num); }
    inline float OpponentDistanceToBall(Unum num) { return PlayerDistanceToBall(TheirSide, num); }

    /* Which side of the field is the object on? */
    inline Fside BallLocationSide() { return LocationSide(BallAbsolutePosition()); }
    inline Fside TeammateLocationSide(Unum t) { return LocationSide(TeammateAbsolutePosition(t)); }
    inline Fside OpponentLocationSide(Unum o) { return LocationSide(OpponentAbsolutePosition(o)); }
    inline Fside PlayerLocationSide(char s, Unum p) { return LocationSide(PlayerAbsolutePosition(s, p)); }

    inline Unum ClosestTeammate()
    {
        if (!MyConf())
            my_error("I'm lost");
        return ClosestTeammateTo(MyPos(), FALSE);
    }
    inline Unum ClosestOpponent()
    {
        if (!MyConf())
            my_error("I'm lost");
        return ClosestOpponentTo(MyPos());
    }

    Unum FurthestBackTeammate(Bool IncludeUs = TRUE, Bool IncludeGoalie = TRUE);
    Unum FurthestBackOpponent();
    /* if we can't find any players, this will return (0,0) */
    Vector PositionOfFurthestBackPlayer(Bool IncludeUs = TRUE);

    Unum FurthestForwardTeammate(Bool IncludeUs = TRUE);
    Unum FurthestForwardOpponent(Bool IncludeGoalie = TRUE);
    /* if we can't find any players, this will return (0,0) */
    Vector PositionOfFurthestForwardPlayer(Bool IncludeUs = TRUE);

    float AngleBetweenClosestTwoOpponents(Vector p);

    Bool InOwnPenaltyArea(Vector p);
    Bool InOwnPenaltyArea();
    Bool BallInOwnPenaltyArea();
    Bool FacingBackInOwnPA();
    Bool FacingBackNearOwnGoal();

    Bool IsPointInBounds(float x, float y, float buffer = 0);
    Bool IsPointInBounds(Vector pt, float buffer = 0) { return IsPointInBounds(pt.x, pt.y, buffer); }
    Rectangle ShiftRectangleInBounds(Rectangle *rect);

    void update();
    void update_self_seen(Time time);
    void update_self(Time time);
    void update_ball(Time time);
    void update_players(Time time);
    void reconcile_unknown_players();

    void VerifyDash(float *dash_power);

    Vector PositionToKickoffPosition(const Vector pos);

    /* offside stuff */

    float my_offside_line;    /* In front of which I'm offside      */
    float their_offside_line; /* Behind which other team is offside */
    void update_offside_lines();

    Bool OffsidePosition(float x, char side);
    inline Bool OffsidePosition(Vector pos, char side) { return OffsidePosition(pos.x, side); }
    Bool TeammateInOffsidePosition(Unum num);
    Bool OpponentInOffsidePosition(Unum num);
    Bool PlayerInOffsidePosition(char side, Unum num);
    inline Bool InOffsidePosition() { return TeammateInOffsidePosition(MyNumber); }
    Unum TeammateOffsideIfIKick();

    float XToAdjustForOffsideX(float x, float buffer);
    inline float XToAdjustForOffsideX(float x) { return XToAdjustForOffsideX(x, CP_at_point_buffer); }
    inline Vector PositionToAdjustForOffsidePosition(Vector pos, float buffer)
    {
        pos.x = XToAdjustForOffsideX(pos.x, buffer);
        return pos;
    }
    inline Vector PositionToAdjustForOffsidePosition(Vector pos)
    {
        pos.x = XToAdjustForOffsideX(pos.x, CP_at_point_buffer);
        return pos;
    }
    Rectangle RectangleToAdjustForOffsideRectangle(Rectangle *rect, float buffer);
    inline Rectangle RectangleToAdjustForOffsideRectangle(Rectangle *rect)
    {
        return RectangleToAdjustForOffsideRectangle(rect, CP_at_point_buffer);
    }

    inline float XToOnsideX(float x, float buffer) { return Min(x, my_offside_line - buffer); }
    inline float XToOnsideX(float x) { return XToOnsideX(x, CP_at_point_buffer); }

    inline Vector PositionToOnsidePosition(Vector pos, float buffer)
    {
        pos.x = XToOnsideX(pos.x, buffer);
        return pos;
    }
    inline Vector PositionToOnsidePosition(Vector pos)
    {
        return PositionToOnsidePosition(pos, CP_at_point_buffer);
    }

    Vector PositionToPullOffsidePosition(Vector pos, float buffer);
    inline Vector PositionToPullOffsidePosition(Vector pos)
    {
        return PositionToPullOffsidePosition(pos, CP_at_point_buffer);
    }

    Bool PullOffside;
    float PullOffsidePosition;

    /* congestion stuff */

    float Congestion(Vector pos, Bool consider_me = FALSE);
    float TeammateCongestion(Unum teammate, Bool consider_me = TRUE);
    inline float MyCongestion() { return TeammateCongestion(MyNumber, FALSE); }
    Unum LeastCongestedTeammate();
    Vector LeastCongestedValidPointInRectangle(Rectangle *rect, Bool attract = FALSE, Vector attract_point = 0);
    Vector LeastCongestedValidPointForPassFromInRectangle(Rectangle *rect, Vector from,
                                                          Bool attract = FALSE, Vector attract_point = 0);

    inline void SetTeammateTired(Unum num, Time time) { TiredTimes[num] = time; }
    inline Bool TeammateTired(Unum num)
    {
        if (num == MyNumber)
            return Tired();
        else
            return (CurrentTime.t - CP_say_tired_interval > SecondLastStartClockTime.t &&
                    CurrentTime - CP_say_tired_interval < TiredTimes[num])
                       ? TRUE
                       : FALSE;
    }

    inline void SetAnnouncedTired(Time time) { SetTeammateTired(MyNumber, time); }
    inline Bool NeedToAnnounceTired()
    {
        if (!Tired())
            my_error("I'm not tired");
        return (CurrentTime.t - CP_say_tired_interval > SecondLastStartClockTime.t &&
                CurrentTime - CP_say_tired_interval < TiredTimes[MyNumber])
                   ? FALSE
                   : TRUE;
    }

    Rectangle OwnPenaltyArea;
    Rectangle OwnGoalieArea;
    Rectangle TheirPenaltyArea;
    Rectangle TheirGoalieArea;
    Rectangle FieldRectangle;
    Vector MyLeftGoalKickSpot;
    Vector MyRightGoalKickSpot;
    Vector TheirLeftGoalKickSpot;
    Vector TheirRightGoalKickSpot;

    /* for the new position based velocity correction */
    Time sight_position_correction_time;
    Vector sight_position_correction;

    /* for ball vcelocity invalidation */
    float quantize_err_const;
    float Tan_of_half_deg; /* the tangent of 1/2 degree */

private:
    Time TiredTimes[MAX_PLAYERS];

    PlayerObject *ClosestPlayerObjectTo(Vector gpos);
    PlayerObject *ClosestTeammateObjectTo(Vector gpos);
    PlayerObject *ClosestOpponentObjectTo(Vector gpos);
    PlayerObject *ClosestPlayerObjectTo(char side, Vector gpos);

    StationaryObject *Fieldline;
    StationaryObject *Marker;
    BallObject Ball;
    PlayerObject *Player;

    SideLine SeenLine;
    MarkerType SeenMarker[MAX_MARKERS];        /* for clearing the seen flag            */
    PlayerObject *MyPlayer[MAX_PLAYERS];       /* pointers to the players on my team    */
    PlayerObject *TheirPlayer[MAX_PLAYERS];    /* pointers to the players on their team */
    PlayerObject *TeamlessPlayer[MAX_PLAYERS]; /* pointers to the players with no team  */
    PlayerObject *FreePlayer[MAX_PLAYERS];     /* elements in Player not assigned       */

    TempPlayerObject *UnknownPlayer; /* Players seen with missing team or num */

    int num_seen_markers;
    int num_my_players;
    int num_their_players;
    int num_teamless_players;
    int num_free_players;
    int num_unknown_players;
    int num_players; /* The maximum -- doesn't include me */
};

/* -*- Mode: C++ -*- */

/* MemAction.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

/****************************************************************************************/

typedef struct PLAYERINTERCEPTINFO
{
    Time time;
    float dash_pow;
    int lookahead;
    InterceptRes res;
    int numCyc;
    Vector pos;
    float dash_pow_to_use;
} PlayerInterceptInfo;

typedef struct TURNKICKCOMMAND
{
    CMDType type;
    Time time;
    float angle;
    float power;
    Bool turn_neck;
    float turn_neck_angle;
} TurnKickCommand;

const int LA_Default = -2;
const int LA_BestSoFar = -1;

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

class ActionInfo : public PositionInfo
{
public:
    void Initialize();

    // these are for hard kicking
    Time HKTime;
    int HKStep;
    int HKStepNext;
    TurnDir HKrot;

    float VelAtPt2VelAtFoot(Vector pt, float targ_vel_at_pt);

    /* these are for ball interception */
    InterceptRes PlayerInterceptionResult(char side, Unum num, float dash_pow);
    InterceptRes PlayerInterceptionResult(char side, Unum num)
    {
        return PlayerInterceptionResult(side, num,
                                        (side == MySide && TeammateTired(num)) ? SP_stamina_inc : SP_max_power);
    }
    Bool PlayerInterceptionAble(char side, Unum num, float dash_pow);
    Bool PlayerInterceptionAble(char side, Unum num)
    {
        return PlayerInterceptionAble(side, num,
                                      (side == MySide && TeammateTired(num)) ? SP_stamina_inc : SP_max_power);
    }
    int PlayerInterceptionNumberCycles(char side, Unum num, float dash_pow);
    int PlayerInterceptionNumberCycles(char side, Unum num)
    {
        return PlayerInterceptionNumberCycles(side, num,
                                              (side == MySide && TeammateTired(num)) ? SP_stamina_inc : SP_max_power);
    }
    Vector PlayerInterceptionPoint(char side, Unum num, float dash_pow);
    Vector PlayerInterceptionPoint(char side, Unum num)
    {
        return PlayerInterceptionPoint(side, num,
                                       (side == MySide && TeammateTired(num)) ? SP_stamina_inc : SP_max_power);
    }
    float PlayerInterceptionDashPower(char side, Unum num, float dash_pow);
    float PlayerInterceptionDashPower(char side, Unum num)
    {
        return PlayerInterceptionDashPower(side, num,
                                           (side == MySide && TeammateTired(num)) ? SP_stamina_inc : SP_max_power);
    }

    InterceptRes TeammateInterceptionResult(Unum num, float dash_pow)
    {
        return PlayerInterceptionResult(MySide, num, dash_pow);
    }
    InterceptRes TeammateInterceptionResult(Unum num)
    {
        return (num == MyNumber) ? MyInterceptionResult() : TeammateInterceptionResult(num, TeammateTired(num) ? SP_stamina_inc : SP_max_power);
    }
    Bool TeammateInterceptionAble(Unum num, float dash_pow)
    {
        return PlayerInterceptionAble(MySide, num, dash_pow);
    }
    Bool TeammateInterceptionAble(Unum num)
    {
        return (num == MyNumber) ? MyInterceptionAble() : TeammateInterceptionAble(num, TeammateTired(num) ? SP_stamina_inc : SP_max_power);
    }
    int TeammateInterceptionNumberCycles(Unum num, float dash_pow)
    {
        return PlayerInterceptionNumberCycles(MySide, num, dash_pow);
    }
    int TeammateInterceptionNumberCycles(Unum num)
    {
        return (num == MyNumber) ? MyInterceptionNumberCycles() : TeammateInterceptionNumberCycles(num, TeammateTired(num) ? SP_stamina_inc : SP_max_power);
    }
    Vector TeammateInterceptionPoint(Unum num, float dash_pow)
    {
        return PlayerInterceptionPoint(MySide, num, dash_pow);
    }
    Vector TeammateInterceptionPoint(Unum num)
    {
        return (num == MyNumber) ? MyInterceptionPoint() : TeammateInterceptionPoint(num, TeammateTired(num) ? SP_stamina_inc : SP_max_power);
    }
    float TeammateInterceptionDashPower(Unum num, float dash_pow)
    {
        return PlayerInterceptionDashPower(MySide, num, dash_pow);
    }
    float TeammateInterceptionDashPower(Unum num)
    {
        return (num == MyNumber) ? MyInterceptionDashPower() : TeammateInterceptionDashPower(num, TeammateTired(num) ? SP_stamina_inc : SP_max_power);
    }

    InterceptRes OpponentInterceptionResult(Unum num, float dash_pow)
    {
        return PlayerInterceptionResult(TheirSide, num, dash_pow);
    }
    InterceptRes OpponentInterceptionResult(Unum num)
    {
        return OpponentInterceptionResult(num, SP_max_power);
    }
    Bool OpponentInterceptionAble(Unum num, float dash_pow)
    {
        return PlayerInterceptionAble(TheirSide, num, dash_pow);
    }
    Bool OpponentInterceptionAble(Unum num)
    {
        return OpponentInterceptionAble(num, SP_max_power);
    }
    int OpponentInterceptionNumberCycles(Unum num, float dash_pow)
    {
        return PlayerInterceptionNumberCycles(TheirSide, num, dash_pow);
    }
    int OpponentInterceptionNumberCycles(Unum num)
    {
        return OpponentInterceptionNumberCycles(num, SP_max_power);
    }
    Vector OpponentInterceptionPoint(Unum num, float dash_pow)
    {
        return PlayerInterceptionPoint(TheirSide, num, dash_pow);
    }
    Vector OpponentInterceptionPoint(Unum num)
    {
        return OpponentInterceptionPoint(num, SP_max_power);
    }
    float OpponentInterceptionDashPower(Unum num, float dash_pow)
    {
        return PlayerInterceptionDashPower(TheirSide, num, dash_pow);
    }
    float OpponentInterceptionDashPower(Unum num)
    {
        return OpponentInterceptionDashPower(num, SP_max_power);
    }

    InterceptRes MyInterceptionResult(float dash_pow)
    {
        return PlayerInterceptionResult(MySide, MyNumber, dash_pow);
    }
    InterceptRes MyInterceptionResult()
    {
        return MyInterceptionResult(CorrectDashPowerForStamina(SP_max_power));
    }
    Bool MyInterceptionAble(float dash_pow)
    {
        return PlayerInterceptionAble(MySide, MyNumber, dash_pow);
    }
    Bool MyInterceptionAble()
    {
        return MyInterceptionAble(CorrectDashPowerForStamina(SP_max_power));
    }
    int MyInterceptionNumberCycles(float dash_pow)
    {
        return PlayerInterceptionNumberCycles(MySide, MyNumber, dash_pow);
    }
    int MyInterceptionNumberCycles()
    {
        return MyInterceptionNumberCycles(CorrectDashPowerForStamina(SP_max_power));
    }
    Vector MyInterceptionPoint(float dash_pow)
    {
        return PlayerInterceptionPoint(MySide, MyNumber, dash_pow);
    }
    Vector MyInterceptionPoint()
    {
        return MyInterceptionPoint(CorrectDashPowerForStamina(SP_max_power));
    }
    float MyInterceptionDashPower(float dash_pow)
    {
        return PlayerInterceptionDashPower(MySide, MyNumber, dash_pow);
    }
    float MyInterceptionDashPower()
    {
        return MyInterceptionDashPower(CorrectDashPowerForStamina(SP_max_power));
    }

    /* just min of what's been done so far - returns -1 if nothing done */
    int GetInterceptionMinCyc();
    inline void ResetInterceptionMinCyc() { IntMinCycTime -= 1; }
    inline int GetInterceptionLookahead() { return InterceptLookahead; }
    void SetInterceptionLookahead(int newval);

    Bool BallPathInterceptValid();
    Vector BallPathInterceptPoint();
    Bool BallPathInterceptAmIThere(float buffer);
    Bool BallPathInterceptAmIThere()
    {
        return BallPathInterceptAmIThere(CP_at_point_buffer);
    }
    float BallPathInterceptDistance();
    /* careful! if ball is kickable, next func returns 0 */
    int BallPathInterceptCyclesForBall();
    Bool BallPathInterceptCanIGetThere(float max_pow = 100.0);

    KickMode BestKickModeAbs(AngleDeg abs_ang);
    KickMode BestKickMode(AngleDeg rel_ang) /* Angle relative to body */
    {
        return BestKickModeAbs(GetNormalizeAngleDeg(rel_ang + MyBodyAng()));
    }

    int EstimatedCyclesToSteal(Unum opp, Vector ball_pos);
    inline int EstimatedCyclesToSteal(Unum opp, AngleDeg ball_ang) // absolute angle
    {
        return EstimatedCyclesToSteal(opp, MyPos() +
                                               Polar2Vector(CP_opt_ctrl_dist, ball_ang));
    }
    inline int EstimatedCyclesToSteal(Unum opp)
    {
        if (!BallPositionValid())
            my_error("EstimateCyclesToSteal: don;t know where ball is");
        return EstimatedCyclesToSteal(opp, BallAbsolutePosition());
    }

    Bool WillDashHelpKick(Vector pt, float dash_pow);
    Bool WillDashHelpKick(Vector pt)
    {
        return WillDashHelpKick(pt, SP_max_power);
    }

    Bool KickInProgress();
    void StartKick(AngleDeg target_angle, KickMode mode, float target_vel, TurnDir rot = TURN_AVOID);
    void StartShot(AngleDeg target_angle, KickMode mode, TurnDir rot = TURN_AVOID);
    void StartPass(Unum target, float target_vel, TurnDir rot = TURN_AVOID);

    AngleDeg kick_in_progress_abs_angle;
    float kick_in_progress_target_vel;
    KickMode kick_in_progress_mode;
    TurnDir kick_in_progress_rotation;

    Unum team_receiver;
    Unum team_passer;
    Time team_pass_time;

    Unum FastestTeammateToBall();
    Unum FastestOpponentToBall();
    Unum BallPossessor(); /* possessor means can get there quickest */
    char TeamInPossession();

    Unum LastBallPossessor;
    Time LastBallPossessorTime;

private:
    Bool kick_in_progress;
    Time start_kick_time;
    Time kick_in_progress_time;

    PlayerInterceptInfo *TeamIntInfo[MAX_PLAYERS];
    PlayerInterceptInfo *OppIntInfo[MAX_PLAYERS];

    Time Stored_Fastest_Teammate_Time;
    Unum Stored_Fastest_Teammate;
    Time Stored_Fastest_Opponent_Time;
    Unum Stored_Fastest_Opponent;

    int InterceptLookahead; /* can either be a positve number or a LA_constant
                   from above */
    int IntMinCyc;
    Time IntMinCycTime;
    void SetIntMinCyc(int newval);

    void BallIntercept_active(float max_pow_to_use, int max_lookahead,
                              char PlayerSide, Unum PlayerNum,
                              PlayerInterceptInfo *pInfo);
    PlayerInterceptInfo CloseBallInterception(float max_pow, int max_lookahead,
                                              Vector vBallPos, Vector vBallVel);
    PlayerInterceptInfo ActiveCanGetThere(float max_pow, int max_lookahead,
                                          Vector vBallPos, Vector vBallVel,
                                          char side, Unum num,
                                          Vector vPlayerPos, Vector vPlayerVel,
                                          float vPlayerAng, int PlayerAngValid,
                                          bool IsThisMe);

    PlayerInterceptInfo *GetPlayerIntInfo(char side, Unum num);
    PlayerInterceptInfo *VerifyIntInfo(char side, Unum num, float dash_pow);

    bool IsSuccessRes(InterceptRes res)
    {
        return (res == BI_CanChase || res == BI_ReadyToKick);
    }

    Time BPItime;
    Bool BPIvalid;
    Bool BPIable;
    float BPIdist;
    Vector BPIpoint;
    float BPIballcyc;

    void VerifyBPIInfo();
    int GetClosestPointToBallPath(Vector *pvPt, float *pNumCycles,
                                  Vector PlayerPos, Vector BallPos,
                                  Vector BallVel);
    void BallIntercept_passive(float max_pow_to_use,
                               PlayerInterceptInfo *pInfo);
};

inline Bool ActionInfo::BallPathInterceptValid()
{
    VerifyBPIInfo();
    return BPIvalid;
}

/* -*- Mode: C++ -*- */

/* Memory.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

class Memory : public ActionInfo
{
public:
    void Initialize(); // depends on the size of the teams
};

/* extern Memory *const Mem; */ /* it's in client.h */

/* client.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

Memory *Mem;

void turn(AngleDeg ang);
void dash(float power);
void kick(float power, AngleDeg dir);
void goalie_catch(AngleDeg dir);
void move(float x, float y);
inline void move(Vector p) { move(p.x, p.y); }
void disconnect();

void turn_neck(AngleDeg ang);
void change_view(Vqual qual, Vwidth width);
inline void change_view(Vqual qual) { change_view(qual, Mem->ViewWidth); }
inline void change_view(Vwidth width) { change_view(Mem->ViewQuality, width); }

/* -*- Mode: C++ -*- */

/* kick.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

void PatsTest_static();
void PatsTest_turn(void);
void PatsTest_kick(void);
void PatsTest_kick2(void);
void PatsTest_conv(void);

int DoTurnKickCommand(TurnKickCommand kickcom);

/* says whether one striaght kick will get the ball to angle dir,
   distance CP_opt_ctrl_dist for current spot */
int is_straight_kick(float dir, float EndDist, float closeMarg);

TurnKickCommand dokick(float ddir, float ddist, float distFactor = 1.0,
                       Bool *pCanKickToTraj = NULL);

TurnDir RotToAvoidOpponent(float abs_dir);
TurnDir RotClosest(float abs_dir);

// relative angle!
KickToRes turnball_kick(AngleDeg target_dir, TurnDir rotate, Bool StopBall,
                        TurnKickCommand *pCom,
                        float EndDist, float closeMarg, float kickFac);
inline KickToRes turnball_kick(AngleDeg target_dir, TurnDir rotate, Bool StopBall,
                               TurnKickCommand *pCom)
{
    return turnball_kick(target_dir, rotate, StopBall, pCom,
                         Mem->CP_opt_ctrl_dist,
                         Mem->CP_closest_margin, Mem->CP_dokick_factor);
}

// relative angle!
KickToRes TurnballTo(AngleDeg target_dir, TurnDir rotate = TURN_AVOID);

#ifdef NOT_USED
int kick_hard_straight(float dir, int step);

int kick_hard_turning(float dir, TurnDir rot, int step);
#endif

TurnDir KickRotationDirectionAbs(AngleDeg abs_ang, TurnDir rot = TURN_AVOID);
inline TurnDir KickRotationDirection(AngleDeg rel_to_body_ang, TurnDir rot = TURN_AVOID)
{
    return KickRotationDirectionAbs(rel_to_body_ang + Mem->MyBodyAng(), rot);
}

TurnKickCommand kick_hard_moderate(AngleDeg abs_dir, float targ_vel);

// int kick_hard_2step(float abs_dir, TurnDir rot, int step, TurnKickCommand* pCom);

int smart_kick_hard_abs(float abs_dir, KickMode mode, float targ_vel,
                        TurnDir rot = TURN_AVOID);
inline int smart_kick_hard_abs(float abs_dir, KickMode mode,
                               TurnDir rot = TURN_AVOID)
{
    return smart_kick_hard_abs(abs_dir, mode,
                               2 * Mem->SP_ball_speed_max, rot);
}

inline int smart_kick_hard(float rel_to_body_dir, KickMode mode, float targ_vel,
                           TurnDir rot = TURN_AVOID)
{
    return smart_kick_hard_abs(rel_to_body_dir + Mem->MyBodyAng(), mode, targ_vel, rot);
}
inline int smart_kick_hard(float rel_to_body_dir, KickMode mode,
                           TurnDir rot = TURN_AVOID)
{
    return smart_kick_hard_abs(rel_to_body_dir + Mem->MyBodyAng(),
                               mode, 2 * Mem->SP_ball_speed_max, rot);
}

/* takes a point and aims to get ball moving at targ_vel at that point
   uses the quickest mode that will get to that speed */
int smart_pass(Vector pt, float targ_vel_at_pt = 1.0, KickMode mode = KM_Moderate, TurnDir rot = TURN_AVOID);

Bool go_to_static_ball(float kick_ang);
inline Bool go_to_static_ball(Vector pt)
{
    return go_to_static_ball((pt - Mem->BallAbsolutePosition()).dir());
}

/* -*- Mode: C++ -*- */

/* dribble.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

void PatsTest_dribble();
void PatsTest_dribble2();

typedef enum DRIBBLERES
{
    DR_None,
    DR_Error,
    DR_LostBall,
    DR_GotThere,
    DR_Going
} DribbleRes;

typedef enum DRIBBLEMODE
{
    DM_None,
    DM_Strict, /* moves very little if ball is at wrong angle */
    DM_Lazy    /* will pick an intermediate point if switching sides, but will
          keep moving. Not as safe with players around */
} DribbleMode;

TurnKickCommand dribble_dash(Vector vEndPos, float max_pow,
                             AngleDeg drib_ang, DribbleMode mode);
TurnKickCommand dribble_kick(Vector vEndPos, float max_pow,
                             AngleDeg drib_ang, DribbleMode mode);
AngleDeg NormalizeDribbleAngle(AngleDeg ang);

/* drib_ang should be between 90 and -90
   we always dribble in front of us */
DribbleRes DribbleTo(Vector vEndPos, float max_dash_pow, float buffer,
                     AngleDeg drib_ang, DribbleMode mode,
                     Bool IsDodge = FALSE, Vector DodgePoint = 0);

DribbleRes SmartDribbleTo(Vector vEndPos, float max_dash_pow, float buffer);
inline DribbleRes SmartDribbleTo(Vector vEndPos, float max_dash_pow)
{
    return SmartDribbleTo(vEndPos, max_dash_pow, Mem->CP_at_point_buffer);
}
inline DribbleRes SmartDribbleTo(Vector vEndPos)
{
    return SmartDribbleTo(vEndPos, Mem->CP_dribble_dash_pow, Mem->CP_at_point_buffer);
}

DribbleRes SimpleDribbleTo(Vector vEndPos, float max_dash_pow = 75,
                           float buffer = 1.0);

/* -*- Mode: C++ -*- */

/* behave.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

void behave();

ActionQueueRes scan_field_with_body();
void turn_neck_to_relative_angle(AngleDeg ang);
void scan_field_with_neck();

ActionQueueRes face_only_body_to_point(Vector point);
void face_only_neck_to_point(Vector point);
ActionQueueRes face_neck_to_point(Vector point);
ActionQueueRes face_neck_and_body_to_point(Vector point);

ActionQueueRes face_only_body_to_player(char side, Unum num);
void face_only_neck_to_player(char side, Unum num);
ActionQueueRes face_neck_to_player(char side, Unum num);
ActionQueueRes face_neck_and_body_to_player(char side, Unum num);

ActionQueueRes face_only_body_to_opponent(Unum opponent);
void face_only_neck_to_opponent(Unum opponent);
ActionQueueRes face_neck_to_opponent(Unum opponent);
ActionQueueRes face_neck_and_body_to_opponent(Unum opponent);

ActionQueueRes face_only_body_to_teammate(Unum teammate);
void face_only_neck_to_teammate(Unum teammate);
ActionQueueRes face_neck_to_teammate(Unum teammate);
ActionQueueRes face_neck_and_body_to_teammate(Unum teammate);

ActionQueueRes face_only_body_to_ball();
void face_only_neck_to_ball();
ActionQueueRes face_neck_to_ball();
ActionQueueRes face_neck_and_body_to_ball();

void get_ball();
void stop_ball();
void hold_ball();
void pass_ball(Unum teammate, float target_vel = 1.0);
void kick_ball(AngleDeg target_angle, KickMode mode, float targ_vel,
               TurnDir rotation = TURN_NONE);
void kick_ball(Vector point, KickMode mode, float targ_vel,
               TurnDir rotation = TURN_NONE);
void kick_ball(AngleDeg target_angle, KickMode mode = KM_Moderate,
               TurnDir rotation = TURN_NONE);
void kick_ball(Vector point, KickMode mode = KM_Moderate,
               TurnDir rotation = TURN_NONE);

ActionQueueRes go_to_point(Vector p, float buffer = 0, float dash_power = 100, DodgeType dodge = DT_all);

/* -*- Mode: C++ -*- */

/* parse.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

void Parse(char *SensoryInfo);

/* -*- Mode: C++ -*- */

/* test.h
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

void test_scan_with_body();
void test_random_movement_in_rectangle(Rectangle *rect);
void test_1v1();
void test_volley();
void test_go_to_ball(AngleDeg kick_angle);
void test_go_to_ball();
void test_go_to_point(Vector p, float buffer, float dash_power = 100);
void test_face_ball();
void test_random_movement();
void test_straight_to_ball();
void test_run_straight();
void test_turn_and_dash_slow();
void test_print_ball();
void test_print_positions();
void test_turnball();
void test_turnball2();
void test_hard_kick(KickMode km);
void test_intercept();
void test_go_to_static_ball();
void test_pred_cycles_to_point();

void test_log_action();

/* -*- Mode: C++ -*- */

/* netif.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

/*
 *  udp client program.
 */

int wait_message(char *buf, Socket *sock)
{
    if (receive_message(buf, sock) == 1)
    {
        while (receive_message(buf, sock) == 1)
            ;
        return 1;
    }
    else
        for (int i = 0; i < 100; i++)
        {
            if (receive_message(buf, sock) == 1)
                return 1;
            // my_error("sleeping, waiting for message");
            usleep(50000);
        }
    return 0;
}

Socket init_connection(char *host, int port)
{
    struct hostent *host_ent;
    struct in_addr *addr_ptr;
    struct sockaddr_in cli_addr;
    int sockfd, val;
    Socket sock;

    sock.socketfd = -1;

    if ((host_ent = (struct hostent *)gethostbyname(host)) == NULL)
    {
        /* Check if a numeric address */
        if ((int)inet_addr(host) == -1)
        {
            return sock;
        }
    }
    else
    {
        addr_ptr = (struct in_addr *)*host_ent->h_addr_list;
        host = inet_ntoa(*addr_ptr);
    }

    sigset_t sigmask;

    sigemptyset(&sigmask);
    sigaddset(&sigmask, SIGIO);
    sigprocmask(SIG_BLOCK, &sigmask, NULL);

    /*
     *  Open UDP socket.
     */
    if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0)
    {
        return sock; /* Can't open socket. */
    }

    if (fcntl(sockfd, F_SETOWN, getpid()) == -1)
        return sock;

    val = fcntl(sockfd, F_GETFL, 0);

/* was "defined NewsOS || defined IRIX" */
#if 1
    val |= O_NDELAY;
#else
    val |= O_NONBLOCK;
#endif
    val |= FASYNC;
    fcntl(sockfd, F_SETFL, val);

    /*
     *  Bind any local address.
     */
    bzero((char *)&cli_addr, sizeof(cli_addr));
    cli_addr.sin_family = AF_INET;
    cli_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    cli_addr.sin_port = htons(0);

    if (bind(sockfd, (struct sockaddr *)&cli_addr,
             sizeof(cli_addr)) < 0)
    {
        return sock; /* Can't bind local address */
    }

    /*
     *  Fill in the structure with the address of the server.
     */
    sock.socketfd = sockfd;

    bzero((char *)&sock.serv_addr, sizeof(sock.serv_addr));
    sock.serv_addr.sin_family = AF_INET;
    sock.serv_addr.sin_addr.s_addr = inet_addr(host);
    sock.serv_addr.sin_port = htons(port);

    return sock;
}

int send_message(char *buf, Socket *sock)
{
    std::cout << "send message: " << buf << std::endl;
    // struct timeval tv_new;
    int n;

    if (!buf)
        return 0;
    n = strlen(buf) + 1;

    if (Mem->CP_save_log && Mem->Initialized)
        fprintf(Mem->SaveLogFile, ">> %s\n\n", buf);

    /* not needed after server 5.23
    gettimeofday(&tv_new, NULL); // no time zone info;
    Mem->SetRealTimeOfLastSend(tv_new);
    */

    if (sendto(sock->socketfd, buf, n, 0,
               (struct sockaddr *)&sock->serv_addr, sizeof(sock->serv_addr)) != n)
        return (-1);
    return 0;
}

int receive_message(char *buf, Socket *sock)
{
    int n;
    socklen_t servlen;
    struct sockaddr_in serv_addr;

    servlen = sizeof(serv_addr);
    n = recvfrom(sock->socketfd, buf, MAXMESG, 0,
                 (struct sockaddr *)&serv_addr, &servlen);
    // std::cout << buf << std::endl;

    if (n < 0)
        if (n == -1 && errno == EWOULDBLOCK)
        {
            buf[0] = '\0';
            /* buf[0] = '\0' ; by Pat- don't want to erase previous messages */
            /* my_error("Receive would block: %d", sock->socketfd); */
            return 0;
        }
        else
        {
            my_error("Receive error: %d on %d", errno, sock->socketfd);
            perror("");
            fflush(stderr);
            return (-1);
        }
    else
    {
        buf[n] = '\0';

        if (Mem->CP_save_log && Mem->Initialized)
        {
            fprintf(Mem->SaveLogFile, "%s\n", buf);
            if (Mem->SaveLogCounter++ % Mem->CP_save_freq == 0)
            {
                fclose(Mem->SaveLogFile);
                Mem->SaveLogFile = fopen(Mem->SaveLogFileName, "a");
            }
        }
        if (Mem->CP_save_sound_log && Mem->Initialized && buf[1] == 'h')
        { /* Hear message */
            fprintf(Mem->SaveSoundLogFile, "%s\n", buf);
            if (Mem->SaveSoundLogCounter++ % Mem->CP_save_sound_freq == 0)
            {
                fclose(Mem->SaveSoundLogFile);
                Mem->SaveSoundLogFile = fopen(Mem->SaveSoundLogFileName, "a");
            }
        }

        sock->serv_addr.sin_port = serv_addr.sin_port;

        if (n == 0)
        {
            my_error("Received null message");
            return 0;
        }
        else
            return 1;
    }
}

void close_connection(Socket *sock)
{
    close(sock->socketfd);
}

/* -*- Mode: C++ -*- */

/* utils.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

int dump_core(char *msg)
{
    my_stamp;
    fprintf(stderr, "dumping core");
    msg[1000000] = 0;
    my_error("Core didn't dump");
    return 0;
}

#define MY_ERROR_LOG_LEVEL 5

void my_error(char *msg)
{
    fprintf(stderr, "PETER's ERROR (player %d, time %d.%d): %s\n",
            Mem->MyNumber, Mem->CurrentTime.t, Mem->CurrentTime.s, msg);
    Mem->LogAction3(MY_ERROR_LOG_LEVEL, "MyError: %s", msg);
    // msg[1000000]=0; /* to segfault */
    if (Mem->CP_stop_on_error)
    {
        int tmp;
        // scanf("%d",&tmp);
    }
}

void my_error(char *msg, int param)
{
    char outstring[100];
    sprintf(outstring, msg, param);
    my_error(outstring);
}

void my_error(char *msg, int param1, int param2)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2);
    my_error(outstring);
}

void my_error(char *msg, int param1, int param2, int param3)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3);
    my_error(outstring);
}

void my_error(char *msg, int param1, int param2, int param3, int param4)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3, param4);
    my_error(outstring);
}

void my_error(char *msg, int param1, int param2, int param3, int param4, int param5)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3, param4, param5);
    my_error(outstring);
}

void my_error(char *msg, int param1, int param2, int param3, int param4, int param5, int param6)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3, param4, param5, param6);
    my_error(outstring);
}

void my_error(char *msg, int param1, int param2, int param3, int param4, int param5, char c1, int param6)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3, param4, param5, c1, param6);
    my_error(outstring);
}

void my_error(char *msg, float param)
{
    char outstring[100];
    sprintf(outstring, msg, param);
    my_error(outstring);
}

void my_error(char *msg, float param1, float param2)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2);
    my_error(outstring);
}

void my_error(char *msg, float param1, float param2, float param3)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3);
    my_error(outstring);
}

void my_error(char *msg, float param1, float param2, float param3, float param4)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3, param4);
    my_error(outstring);
}

void my_error(char *msg, float param1, int param2)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2);
    my_error(outstring);
}

void my_error(char *msg, char *param)
{
    char outstring[100];
    sprintf(outstring, msg, param);
    my_error(outstring);
}

void my_error(char *msg, char param1, int param2, int param3)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3);
    my_error(outstring);
}

void my_error(char *msg, char param1, int param2, float param3, float param4)
{
    char outstring[100];
    sprintf(outstring, msg, param1, param2, param3, param4);
    my_error(outstring);
}

float int_pow(float x, int p)
{
    if (p < 0)
        return (1.0 / int_pow(x, -p));
    else
    {
        float ans = 1.0;
        for (int i = 0; i < p; i++)
            ans *= x;
        return ans;
    }
}

int closest_int(float x)
{
    if (x < 0)
        x -= 0.5;
    else
        x += 0.5;
    return (int)x;
}

void NormalizeAngleDeg(int *ang)
{
    if (fabs(*ang) > 5000)
        my_error("Huge angle passed to NormalizeAngleDeg");
    while (*ang > 180)
        *ang -= 360;
    while (*ang < -180)
        *ang += 360;
}

void NormalizeAngleDeg(AngleDeg *ang)
{
    if (fabs(*ang) > 5000)
        my_error("Huge angle passed to NormalizeAngleDeg");
    while (*ang > 180)
        *ang -= 360;
    while (*ang < -180)
        *ang += 360;
}

void NormalizeAngleRad(AngleRad *ang)
{
    while (*ang > M_PI)
        *ang -= 2 * M_PI;
    while (*ang < -M_PI)
        *ang += 2 * M_PI;
}

AngleDeg GetNormalizeAngleDeg(AngleDeg ang)
{
    if (fabs(ang) > 5000)
        my_error("Huge angle passed to GetNormalizeAngleDeg");
    while (ang > 180)
        ang -= 360;
    while (ang < -180)
        ang += 360;
    return ang;
}

float GetDistance(float *x, float *y, float *a, float *b)
{
    return sqrt((*x - *a) * (*x - *a) + (*y - *b) * (*y - *b));
}

float weighted_avg(float val1, float val2, float w1, float w2)
{
    return (val1 * w1 + val2 * w2) / (w1 + w2);
}

float SumGeomSeries(float first_term, float r, int n)
{
    return first_term * (Exp(r, n) - 1) / (r - 1);
}

float SolveForLengthGeomSeries(float first_term, float r, float sum)
{
    if (r < 0)
        my_error("SolveForSumGeomSeries: can't take r < 0");
    float temp = sum * (r - 1) / first_term + 1;
    if (temp <= 0)
        return -1.0;
    return log(temp) / log(r);
}

float SolveForFirstTermGeomSeries(float r, int n, float sum)
{
    return sum * (r - 1) / (Exp(r, n) - 1);
}

/* returns a pointer to a static buffer, so be careful! */
char *repeat_char(char c, int n)
{
    const int MAX_REP = 100;
    static char out[MAX_REP + 1];
    if (n > MAX_REP)
        my_error("repeat_char: asking for too many characters");
    for (int i = 0; i < n; i++)
        out[i] = c;
    out[n] = 0;
    return out;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

Time Time::operator-(const int &a)
{
    if (s == 0 && t - a > Mem->LastStartClockTime.t) /* default case */
        return Time(t - a, 0);

    if (s > 0)
    {
        if (a <= s)
            return Time(t, s - a); /* Just take off from stopped time */
        else
        {                                         /* a > s */
            Time new_time = Time(t - (a - s), 0); /* take off from stopped time, then server time */
            if (new_time < Mem->LastStartClockTime)
                my_error("can't be sure of this subtraction (1)");
            return new_time;
        }
    }
    else
    { /* t-a <= Mem->LastStartClockTime.t */
        int stopped_time = a - (t - Mem->LastStartClockTime.t);
        if (stopped_time > Mem->LastStartClockTime.s)
        {
            if (!Mem->LastStartClockTime.t) /* if ==0, OK---account for players starting at different times */
                return Time(0, 0);
            int tmp = Mem->LastStartClockTime.t - (stopped_time - Mem->LastStartClockTime.s);
            if (tmp <= Mem->SecondLastStartClockTime.t)
                my_error("can't be sure of this subtraction (2) %d %d %d  %d", t, s, a, Mem->LastStartClockTime.s);
            return Time(tmp, 0);
        }
        return Time(Mem->LastStartClockTime.t, Mem->LastStartClockTime.s - stopped_time);
    }
}

int Time::operator-(const Time &a)
{
    if (s == 0)
    {
        if (a.t < Mem->LastStartClockTime.t)
        {
            if (a.t <= Mem->SecondLastStartClockTime.t)
                my_error("Can't be sure of this subtraction (3) %d %d  -  %d %d", t, s, a.t, a.s);
            return (t - a.t) + Mem->LastStartClockTime.s;
        }
        if (a.t > Mem->LastStartClockTime.t)
            return (t - a.t);
        if (a.t == Mem->LastStartClockTime.t)
            return (t - a.t) + (Mem->LastStartClockTime.s - a.s);
    }
    else if (s > 0)
    {
        if (a.t <= Mem->SecondLastStartClockTime.t) /* If they're equal, it's not OK */
            my_error("Can't be sure of this subtraction (4) %d %d  %d %d", t, s, a.t, a.s);
        if (a.t < Mem->LastStartClockTime.t)
            return Mem->LastStartClockTime.s + (t - a.t);
        else if (a.t == Mem->LastStartClockTime.t && a.t != t) /* a is during the last stopped interval */
            return (s + (Mem->LastStartClockTime.s - a.s)) + (t - a.t);
        return (s - a.s) + (t - a.t);
    }
    else /* s<0 */
        my_error("s shouldn't be <0");
    return 0;
}

Time Time::operator+(const int &a)
{
    if (s == 0 && t > Mem->LastStartClockTime.t && t + a < Mem->CurrentTime.t) /* default case */
        return Time(t + a, 0);

    Time new_time;

    if (s == 0)
    {
        if (Mem->LastStartClockTime.t > t)
        { /* Could've missed one already */
            my_error("Can't be sure of this addition (1) %d %d", Mem->LastStartClockTime.t, t);
            new_time = Time(t + a, 0);
        }
        if (t + a < Mem->CurrentTime.t)
            new_time = Time(t + a, 0);
        else /* t+a >= Mem->CurrentTime.t */
            new_time = Time(Mem->CurrentTime.t, a - (Mem->CurrentTime.t - t));
    }
    else if (s > 0)
    {
        if (t == Mem->CurrentTime.t) /* clock still stopped */
            new_time = Time(t, s + a);
        else
        {
            if (Mem->LastStartClockTime.t != t) /* Could've missed one already */
                my_error("Can't be sure of this addition (2)");
            new_time = Time(t + (a - (Mem->LastStartClockTime.s - s)), 0);
        }
    }
    else /* s<0 */
        my_error("s shouldn't be <0");

    if (new_time > Mem->CurrentTime) /* clock might stop */
        my_error("Can't be sure of this addition (3)");

    return new_time;
}

Bool Time::CanISubtract(const Time &a)
{
    if (s == 0)
    {
        if (a.t < Mem->LastStartClockTime.t)
        {
            if (a.t <= Mem->SecondLastStartClockTime.t)
                return FALSE;
            return TRUE;
        }
        return TRUE;
    }
    else if (s > 0)
    {
        if (a.t <= Mem->SecondLastStartClockTime.t) /* If they're equal, it's not OK */
            return FALSE;
        return TRUE;
    }
    else /* s<0 */
        my_error("s shouldn't be <0");
    return FALSE;
}

/****************************************************************************/
/* These routines are to save time instead of using sscanf or atof, etc.    */
/* When passing **str_ptr, the pointer is advanced past the number          */
/* When passing  *str    , the pointer remains where it was before          */
/****************************************************************************/

double get_double(char **str_ptr)
{

    double d_frac, result;
    int m_flag, d_flag;

    d_frac = 1.0;
    result = 0.0;
    m_flag = d_flag = 0;

    /* Advance to the beginning of the number */
    while (!isdigit(**str_ptr) && **str_ptr != '-' && **str_ptr != '.')
        (*str_ptr)++;

    /* Process the number bit by bit */
    while ((**str_ptr != ')') && (**str_ptr != ' ') && (**str_ptr != '\n') && (**str_ptr != ','))
    {
        if (**str_ptr == '.')
            d_flag = 1;
        else if (**str_ptr == '-')
            m_flag = 1;
        else if (d_flag)
        {
            d_frac *= 10.0;
            result += (double)(**str_ptr - '0') / d_frac;
        }
        else
            result = result * 10.0 + (double)(**str_ptr - '0');
        (*str_ptr)++;
    }
    if (m_flag)
        result = -result;

    // printf("%.1f\n",result);

    return result;
}

/* Get the number, but don't advance pointer */

double get_double(char *str)
{
    char **str_ptr = &str;
    return get_double(str_ptr);
}

/****************************************************************************/

float get_float(char **str_ptr)
{

    float d_frac, result;
    int m_flag, d_flag;

    d_frac = 1.0;
    result = 0.0;
    m_flag = d_flag = 0;

    /* Advance to the beginning of the number */
    while (!isdigit(**str_ptr) && **str_ptr != '-' && **str_ptr != '.')
        (*str_ptr)++;

    /* Process the number bit by bit */
    while ((**str_ptr != ')') && (**str_ptr != ' ') && (**str_ptr != '\n') && (**str_ptr != ','))
    {
        if (**str_ptr == 'e')
            my_error("There's an e in my float! %s", *str_ptr);
        if (**str_ptr == '.')
            d_flag = 1;
        else if (**str_ptr == '-')
            m_flag = 1;
        else if (d_flag)
        {
            d_frac *= 10.0;
            result += (float)(**str_ptr - '0') / d_frac;
        }
        else
            result = result * 10.0 + (float)(**str_ptr - '0');
        (*str_ptr)++;
    }
    if (m_flag)
        result = -result;

    // printf("%.1f\n",result);

    return result;
}

/* Get the number, but don't advance pointer */

float get_float(char *str)
{
    char **str_ptr = &str;
    return get_float(str_ptr);
}

/****************************************************************************/

int get_int(char **str_ptr)
{

    int result;
    int m_flag;

    result = 0;
    m_flag = 0;

    /* Advance to the beginning of the number */
    while (!isdigit(**str_ptr) && **str_ptr != '-')
        (*str_ptr)++;

    /* Process the number bit by bit */
    while ((**str_ptr != ')') && (**str_ptr != ' ') && (**str_ptr != '\n') && (**str_ptr != ','))
    {
        if (**str_ptr == '-')
            m_flag = 1;
        else
            result = result * 10 + (int)(**str_ptr - '0');
        (*str_ptr)++;
    }
    if (m_flag)
        result = -result;

    return result;
}

int get_int(char *str)
{
    char **str_ptr = &str;
    return get_int(str_ptr);
}

/****************************************************************************/

void get_word(char **str_ptr)
{
    while (!isalpha(**str_ptr))
        (*str_ptr)++;
}

/****************************************************************************/

void get_next_word(char **str_ptr)
{
    while (isalpha(**str_ptr))
        (*str_ptr)++;
    get_word(str_ptr);
}

/****************************************************************************/

void get_token(char **str_ptr)
{
    advance_past_space(str_ptr);
    while ((*str_ptr) && !isspace(**str_ptr))
        (*str_ptr)++;
}

/****************************************************************************/

void advance_to(char c, char **str_ptr)
{
    while (**str_ptr != c)
        (*str_ptr)++;
}

/****************************************************************************/

void advance_past_space(char **str_ptr)
{
    while ((*str_ptr) && isspace(**str_ptr))
        (*str_ptr)++;
}

/****************************************************************************/
/* These routines are to save time instead of using sprintf or atof, etc.   */
/* *str should point to the END of the string where the number is going     */
/* return the length of the number placed in                                */
/****************************************************************************/

int put_float(char *str, float fnum, int precision)
{
    int m_flag = 0, length = 0;
    int num, old_num;

    for (int i = 0; i < precision; i++)
        fnum *= 10;

    num = closest_int(fnum); /* round off the rest */

    if (precision == 0)
        return put_int(str, num);

    if (num < 0)
    {
        m_flag = 1;
        num = -num;
    }

    old_num = num;
    while (num > 0 || length < precision)
    {
        num /= 10;
        *str = '0' + old_num - num * 10;
        old_num = num;
        str--;
        length++;
        if (length == precision)
        {
            *str = '.';
            str--;
            length++;
            if (num == 0)
            {
                *str = '0';
                str--;
                length++;
                break;
            }
        }
    }

    if (m_flag)
    {
        *str = '-';
        length++;
    }

    return length;
}

/****************************************************************************/

int put_int(char *str, int num)
{

    int m_flag = 0, length = 0;
    int old_num;

    if (num == 0)
    {
        *str = '0';
        return 1;
    }

    if (num < 0)
    {
        m_flag = 1;
        num = -num;
    }

    old_num = num;
    while (num > 0)
    {
        num /= 10;
        *str = '0' + old_num - num * 10;
        old_num = num;
        str--;
        length++;
    }

    if (m_flag)
    {
        *str = '-';
        length++;
    }

    return length;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void BubbleSort(int length, int *elements, float *keys)
{

    /*  for (int a=0; a<length; a++){
        printf("%d--%.1f ",elements[a], keys[a]);
      }
      printf("\n");*/

    /* Sort the elements in increasing order according to the keys */
    float keytemp;
    int eltemp;
    for (int i = 0; i < length; i++)
    {
        for (int j = i + 1; j < length; j++)
        {
            if (keys[j] < keys[i])
            {
                keytemp = keys[i];
                keys[i] = keys[j];
                keys[j] = keytemp;
                eltemp = elements[i];
                elements[i] = elements[j];
                elements[j] = eltemp;
            }
        }
    }
    /*  for (int a=0; a<length; a++){
        printf("%d--%.1f ",elements[a], keys[a]);
      }
      printf("\n");*/
}

/****************************************************************************/

int BinarySearch(int length, float *elements, float key)
{

    /* Assume the list is already sorted in increasing order */
    int lbound = 0, ubound = length;

    for (int index = length / 2; ubound - lbound > 0; index = lbound + (ubound - lbound) / 2)
    {
        /* printf("%d ",index); */
        if (elements[index] == key)
        {
            lbound = ubound = index;
        }
        else if (elements[index] < key)
        {
            lbound = index + 1;
        }
        else
        {
            ubound = index - 1;
        }
    }

    int toReturn = Max(ubound, lbound);
    if (elements[toReturn] < key)
        toReturn++; /* Guarantees >= key */

    return toReturn;
}

/****************************************************************************/

/* replace all occurrences in a string */
void StrReplace(char *str, char oldchar, char newchar)
{
    int i = 0;
#if 0
  int numReplaced;
#endif
    int strLength = strlen(str);
    while (i++ < strLength)
    {
        if (str[i] == oldchar)
        {
            str[i] = newchar;
#if 0
      numReplaced++;
#endif
        }
        if (i == 1000)
            my_error("String of length >1000?");
    }
#if 0
  printf("***Replaced %d %c's in string of length %d (%d): %s***\n",
	 numReplaced,oldchar,strlen(str),i,str);
#endif
}

/****************************************************************************/
/***************************   random stuff    ******************************/
/****************************************************************************/
/* From Andrew's C package                                                  */

int int_random(int n)
{
    static int FirstTime = TRUE;

    if (FirstTime)
    {
        /* initialize the random number seed. */
        timeval tp;
        gettimeofday(&tp, NULL);
        srandom((unsigned int)tp.tv_usec);
        FirstTime = FALSE;
    }

    if (n > 2)
        return (random() % n);
    else if (n == 2)
        return (((random() % 112) >= 56) ? 0 : 1);
    else if (n == 1)
        return (0);
    else
    {
        printf("int_random(%d) ?\n", n);
        my_error("You called int_random(<=0)");
        return (0);
    }
}

float range_random(float lo, float hi)
{
    int x1 = int_random(10000);
    int x2 = int_random(10000);
    float r = (((float)x1) + 10000.0 * ((float)x2)) / (10000.0 * 10000.0);
    return (lo + (hi - lo) * r);
}

int very_random_int(int n)
{
    int result = (int)range_random(0.0, (float)n); /* rounds down */
    if (result == n)
        result = n - 1;
    return (result);
}

void GetStampedName(char *name, char *outputName)
{
    char date[100], weekday[10], month[10], temp[10];
    int day, hour, min, sec, year;
    FILE *dateFile;

    // system("date > ../date.log");        /* Put the date in a file */ /* done by a player */

    dateFile = fopen("../date.log", "r");
    fscanf(dateFile, "%[^\n]", date); /* Scan it in             */
    fclose(dateFile);

    sscanf(date, "%s %s %d %d:%d:%d %s %d",
           weekday, month, &day, &hour, &min, &sec, temp, &year);
    sprintf(name, "%s-%s%d-%d:%d:%d.dat", outputName, month, day, hour, min, sec);
}

/* -*- Mode: C++ -*- */

/* geometry.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

#ifdef DEBUG_OUTPUT
#define DebugGeom(x)
#else
#define DebugGeom(x)
#endif

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        Rectangle Class                                       */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

Rectangle::Rectangle()
{
    left_x = right_x = top_y = bottom_y = 0;
}

/********************************************************************************/

Rectangle::Rectangle(const float l, const float r, const float t, const float b)
{
    if (l >= r)
        my_error("width must be positive");
    if (t >= b)
        my_error("height must be positive");

    left_x = l;
    right_x = r;
    top_y = t;
    bottom_y = b;
}

/********************************************************************************/

Rectangle::Rectangle(const Vector center, const Vector size)
{
    left_x = center.x - size.x / 2.0;
    right_x = center.x + size.x / 2.0;
    top_y = center.y - size.y / 2.0;
    bottom_y = center.y + size.y / 2.0;
}

/********************************************************************************/

Bool Rectangle::IsWithin(const Vector &p)
{
    return (Bool)((p.x >= left_x) && (p.x <= right_x) && (p.y >= top_y) && (p.y <= bottom_y));
}

/********************************************************************************/

Vector Rectangle::AdjustToWithin(const Vector &p)
{

    Vector r = p;

    if (r.y - top_y < 0)
        r.y = top_y;
    if (bottom_y - r.y < 0)
        r.y = bottom_y;
    if (right_x - r.x < 0)
        r.x = right_x;
    if (r.x - left_x < 0)
        r.x = left_x;

    return r;
}

/********************************************************************************/
// order: top, right, bot, left
Line Rectangle::GetEdge(int n)
{
    switch (n % 4)
    {
    case 0:
        return TopEdge();
    case 1:
        return RightEdge();
    case 2:
        return BottomEdge();
    case 3:
        return LeftEdge();
    }
    my_error("Rectangle::GetEdge: how did I get here");
    return Line();
}

/********************************************************************************/
// order: TL, TR, BR, BL
Vector Rectangle::GetPoint(int n)
{
    switch (n % 4)
    {
    case 0:
        return TopLeftCorner();
    case 1:
        return TopRightCorner();
    case 2:
        return BottomRightCorner();
    case 3:
        return BottomLeftCorner();
    }
    my_error("Rectangle::GetPoint: how did I get here");
    return Vector(0, 0);
}

/********************************************************************************/

Vector Rectangle::nearestHEdge(const Vector &p)
/* find nearest horizontal line */
{
    static Vector r;

    r.x = Min(Max(p.x, left_x), right_x);
    r.y = ((p.y - top_y) < (bottom_y - p.y)) ? top_y : bottom_y;

    return r;
}

/********************************************************************************/

Vector Rectangle::nearestVEdge(const Vector &p)
/* find nearest vertical line */
{
    static Vector r;

    r.x = ((p.x - left_x) < (right_x - p.x)) ? left_x : right_x;
    r.y = Min(Max(p.y, top_y), bottom_y);

    return r;
}

/********************************************************************************/

Vector Rectangle::nearestEdge(const Vector &p)
{
    if (Min((p.x - left_x), (right_x - p.x)) < Min((p.y - top_y), (bottom_y - p.y)))
        return nearestVEdge(p);
    else
        return nearestHEdge(p);
}

/********************************************************************************/

Line Rectangle::nearestHEdgeLine(const Vector &p)
{
    return ((p.y - top_y) < (bottom_y - p.y)) ? TopEdge() : BottomEdge();
}

/********************************************************************************/

Line Rectangle::nearestVEdgeLine(const Vector &p)
{
    return ((p.x - left_x) < (right_x - p.x)) ? LeftEdge() : RightEdge();
}

/********************************************************************************/

Line Rectangle::nearestEdgeLine(const Vector &p)
{
    if (Min((p.x - left_x), (right_x - p.x)) < Min((p.y - top_y), (bottom_y - p.y)))
        return nearestVEdgeLine(p);
    else
        return nearestHEdgeLine(p);
}

/********************************************************************************/

float Rectangle::DistanceToEdge(const Vector &p)
{
    if (!IsWithin(p))
    {
        Vector q = AdjustToWithin(p);
        return -(q.dist(p)); /* distance outside is a negative number */
    }

    return Min((p.x - left_x), Min((right_x - p.x), Min((p.y - top_y), (bottom_y - p.y))));
}

/********************************************************************************/

Vector Rectangle::random()
{
    static Vector r;

    r.x = range_random(left_x, right_x);
    r.y = range_random(bottom_y, top_y);

    return r;
}

/********************************************************************************/

Rectangle Rectangle::expand(float val)
{
    return Rectangle(left_x - val, right_x + val, top_y - val, bottom_y + val);
}

/********************************************************************************/

void Rectangle::Print()
{
    printf("RECTANGLE:  x = %.1f to %.1f   y = %.1f to %.1f\n",
           LeftX(), RightX(), TopY(), BottomY());
}

/********************************************************************************/

Vector Rectangle::RayIntersection(Ray r)
{
    if (!IsWithin(r.origin))
    {
        my_error("Rectangle/Ray intersection only handles rays inside the rectangle!");
        return 0;
    }

    Vector int_pt[2];
    int num_int = 0;
    for (int i = 0; i < 4; i++)
    {
        Vector pt;
        if (GetEdge(i).RayIntersection(r, &pt))
            int_pt[num_int++] = pt;
    }

    if (num_int == 0)
    {
        my_error("Rectangle ray intersection: no edge intersection?");
        return 0;
    }
    else if (num_int == 1)
    {
        return int_pt[0]; // same slope as one pair of edges
    }
    else if (num_int == 2)
    {
        Rectangle temp_rect = expand(FLOAT_EPS);
        for (int j = 0; j < 2; j++)
            if (temp_rect.IsWithin(int_pt[j]))
                return int_pt[j];
        my_error("Rectangle ray intersection: how did I get here");
        return 0;
    }
    else
    {
        my_error("Rectangle ray intersection: how did num_int get so big: %d", num_int);
        return 0;
    }
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                         Line Class                                           */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

Line::Line(float x_coef, float y_coef, float constant)
{
    A = x_coef;
    B = y_coef;
    C = constant;
}

/********************************************************************************/

void Line::LineFromTwoPoints(Vector pt1, Vector pt2)
{
    float temp = (pt2.x - pt1.x);
    if (fabs(temp) < FLOAT_EPS)
    {
        if (fabs(pt2.y - pt1.y) < FLOAT_EPS)
            my_error("LineFromTwoPoints: points can not be the same!");
        A = 1;
        B = 0;
    }
    else
    {
        float m = (pt2.y - pt1.y) / temp;
        A = -m;
        B = 1;
    }
    C = -(A * pt2.x + B * pt2.y);
}

/********************************************************************************/

void Line::LineFromRay(Ray r)
{
    if (fabs(r.direction.y) < FLOAT_EPS && fabs(r.direction.x) < FLOAT_EPS)
        my_error("LineFromRay: dir can not be zero");
    LineFromTwoPoints(r.origin, r.origin + r.direction);
}

/********************************************************************************/

Bool Line::PointOnLine(float x, float y)
{
    return (Bool)(fabs(A * x + B * y + C) < FLOAT_EPS);
}

/********************************************************************************/

float Line::dist(Vector pt)
{
    return fabs((A * pt.x + B * pt.y + C) / sqrt(Sqr(A) + Sqr(B)));
}

/********************************************************************************/

float Line::dist2(Vector pt)
{
    return fabs(Sqr(A * pt.x + B * pt.y + C) / (Sqr(A) + Sqr(B)));
}

/********************************************************************************/

float Line::angle()
{
    return ATan(-A / B);
}

/********************************************************************************/

Vector Line::ProjectPointUsingCircle(Vector pt)
{
    /* here's the idea:
       first get the plane equation for the ray
       then use the closest distance formula to get the closest distance
       then using the eq for that line and the eq for a circle of the dist radius
         around our curretn position, find the nearest point.
         Whew! */

    /* compute the dist */
    Vector retPt;
    float d = dist(pt); /* the min distance */

    /* intersect the circle and the line */
    float a, b, c; /* coefficent in quadratic to solve for x */
    float disc;    /* discriminant in quadratic equation */
    a = 1 + Sqr(A);
    b = 2 * (A * (pt.y + C) - pt.x);
    c = Sqr(pt.x) + Sqr(pt.y + C) - Sqr(d);
    disc = Sqr(b) - 4 * a * c;
    /* the discriminant should be zero since this is the radius
       is the closest distance between the center and the line */
    if (fabs(disc) > FLOAT_EPS)
        fprintf(stderr, "GetClosestPointToBallPath: discrimannt is bad! %f\n", disc);
    retPt.x = -b / (2 * a);

    /* we compute two possible solutions for y and then see which one is on the
       line of the ray's path */
    float sol1, sol2;
    sol1 = sqrt(Sqr(d) - Sqr(retPt.x - pt.x));
    sol2 = -sol1;
    sol1 += pt.y;
    sol2 += pt.y;

    if (fabs(A * (retPt.x) + B * sol1 + C) < FLOAT_EPS)
    {
        /* sol1 is on line */
        retPt.y = sol1;
    }
    else if (fabs(A * (retPt.x) + B * sol2 + C) < FLOAT_EPS)
    {
        /* sol2 is on line */
        retPt.y = sol2;
    }
    else
        fprintf(stderr, "GetClosestPointToBallPath: neither solution works!\n");

    DebugGeom(printf("  dist: %f\t ptMod: %f\n", d,
                     sqrt(Sqr(pt.x - retPt.x) + Sqr(pt.y - retPt.y))));

    return retPt;
}

/********************************************************************************/

float Line::get_y(float x)
{
    if (B != 0)
        return (-A * x - C) / B;

    my_error("can't get y");
    return 0;
}

/********************************************************************************/

float Line::get_x(float y)
{
    if (A != 0)
        return (-B * y - C) / A;

    my_error("can't get x");
    return 0;
}

/********************************************************************************/

Bool Line::InBetween(Vector pt, Vector end1, Vector end2)
{
    if (!OnLine(end1) || !OnLine(end2))
        my_error("Line::InBetween: passed in points that weren't on line");

    pt = ProjectPoint(pt);
    float dist2 = end1.dist2(end2);

    return (pt.dist2(end1) <= dist2 && pt.dist2(end2) <= dist2) ? TRUE : FALSE;
}

/********************************************************************************/

Vector Line::GetClosestPtInBetween(Vector pt, Vector end1, Vector end2)
{
    if (!OnLine(end1) || !OnLine(end2))
        my_error("Line::InBetween: passed in points that weren't on line");

    if (InBetween(pt, end1, end2))
        return ProjectPoint(pt);
    if (end1.dist2(pt) < end2.dist2(pt))
        return end1;
    else
        return end2;
}

/********************************************************************************/

Vector Line::intersection(Line l)
{
    Vector result = 0;
    if (SameSlope(l))
    {
        // if ( B == 0 && l.B == 0 || A/B == l.A/l.B ){
        my_error("Lines have same slope");
        DebugGeom(cout << "Lines have same slope" << endl);
        return result;
    }

    if (B == 0)
    {
        result.x = -C / A;
        result.y = l.get_y(result.x);
        return result;
    }

    if (l.B == 0)
    {
        result.x = -l.C / l.A;
        result.y = get_y(result.y);
        return result;
    }

    result.x = (C * l.B - B * l.C) / (l.A * B - A * l.B);
    result.y = get_y(result.x);
    return result;
}

/********************************************************************************/
Bool Line::RayIntersection(Ray r, Vector *ppt)
{
    Line lRay(r);

    if (SameSlope(lRay))
        return FALSE;

    *ppt = intersection(lRay);
    return (r.InRightDir(*ppt))
               // fabs(GetNormalizeAngleDeg((*ppt - r.origin).dir() - r.direction.dir())) < 10)
               ? TRUE
               : FALSE;
}

/********************************************************************************/
Bool Line::IsPtCloserToPtOnLine(Vector pt1, Vector pt2, Vector targ_pt)
{
    if (!OnLine(targ_pt))
        my_error("IsPtCloserToPtOnLine: targ_pt not on line");

    pt1 = ProjectPoint(pt1);
    pt2 = ProjectPoint(pt2);

    return (pt1.dist(targ_pt) < pt2.dist(targ_pt)) ? TRUE : FALSE;
}

/********************************************************************************/

// return TRUE on top/left part of plane
Bool Line::HalfPlaneTest(Vector pt)
{
    if (B == 0)
        return (pt.x < -C / A) ? TRUE : FALSE;
    return (pt.y > get_y(pt.x)) ? TRUE : FALSE;
}

Bool Line::SameSlope(Line l)
{
    return (B == 0 && l.B == 0 || A / B == l.A / l.B) ? TRUE : FALSE;
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                       Ray Class                                              */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
Ray::Ray(Vector orig, Vector dir)
{
    origin = orig;
    if (fabs(dir.y) < FLOAT_EPS && fabs(dir.x) < FLOAT_EPS)
    {
        my_error("Ray: dir can not be zero");
        direction = Vector(1, 0);
    }
    else
    {
        direction = dir;
        direction = direction.Normalize();
    }
}

Bool Ray::OnRay(Vector pt)
{
    Vector v = pt - origin;
    return (fabs(Sin(v.dir() - direction.dir()) * v.mod()) < FLOAT_EPS)
               ? TRUE
               : FALSE;
}

Bool Ray::InRightDir(Vector pt)
{
    return (fabs(GetNormalizeAngleDeg((pt - origin).dir() - direction.dir())) < 10)
               ? TRUE
               : FALSE;
}

Bool Ray::intersection(Line l, Vector *pPt)
{
    return l.RayIntersection(*this, pPt);
}

Bool Ray::intersection(Ray r, Vector *pPt)
{
    Line thisLine(*this), argLine(r);

    if (thisLine.SameSlope(argLine))
        return FALSE;
    *pPt = thisLine.intersection(argLine);

    /* now make sure that the intersection is the correct direction on both lines */
    return (InRightDir(*pPt) && r.InRightDir(*pPt))
               //	   fabs(GetNormalizeAngleDeg((*pPt - origin).dir() - direction.dir())) < 10 &&
               // fabs(GetNormalizeAngleDeg((*pPt - r.origin).dir() - r.direction.dir())) < 10)
               ? TRUE
               : FALSE;
}

/* intersects a ray and a cricle */
/* return the number of solutions */
/* psol1 1 is not as afar along the ray as psol2 */
int Ray::CircleIntersect(float rad, Vector center, Vector *psol1, Vector *psol2)
{
    DebugGeom(cout << "RCI: origin: " << origin << "\tdirection: " << direction << endl
                   << "rad: " << rad << "\tcenter: " << center << endl);
    float a, b, c, disc;
    float t1, t2;
    a = Sqr(direction.x) + Sqr(direction.y);
    b = 2.0 * ((origin.x - center.x) * direction.x + (origin.y - center.y) * direction.y);
    c = Sqr(origin.x - center.x) + Sqr(origin.y - center.y) - Sqr(rad);
    DebugGeom(printf(" RCI: a: %f\tb: %f\t c: %f\n", a, b, c));

    disc = Sqr(b) - 4 * a * c;
    if (disc < 0)
    {
        DebugGeom(printf(" RCI disc < 0: %f\n", disc));
        return 0;
    }

    disc = sqrt(disc);
    t1 = (-b + disc) / (2.0 * a);
    t2 = (-b - disc) / (2.0 * a);
    DebugGeom(printf(" RCI: t1: %f\tt2: %f\n", t1, t2));

    if (t1 > t2)
    {
        DebugGeom(printf(" RCI: reversing t1, t2\n"));
        float temp = t1;
        t1 = t2;
        t2 = temp;
    }

    if (t1 > 0.0)
    {
        if (t2 > 0.0)
        {
            *psol1 = origin + direction * t1;
            *psol2 = origin + direction * t2;
            DebugGeom(printf(" RCI:two sols\n"));
            return 2;
        }
        else
        {
            my_error("RayCircleIntersect: weird roots");
            return 0;
        }
    }
    else if (t2 > 0.0)
    {
        *psol1 = origin + direction * t2;
        DebugGeom(printf(" RCI:t2 only sol\n"));
        return 1;
    }
    else
        return 0;

    return 0;
}

Vector Ray::RectangleIntersection(Rectangle R)
{
    return R.RayIntersection(*this);
}

Vector Ray::GetClosestPoint(Vector pt)
{
    Line l(*this);
    Vector close_pt = l.ProjectPoint(pt);
    if (OnRay(close_pt))
        return close_pt;
    else
        return origin;
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                  Miscellaneous Geometry Functions                            */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

#ifdef OLD_CODE
/* intersects a ray and a cricle */
/* return the number of solutions */
/* psol1 1 is not as afar along the ray as psol2 */
int RayCircleIntersect(Ray r, float rad, Vector center,
                       Vector *psol1, Vector *psol2)
{
    DebugGeom(cout << "RCI: r.origin: " << r.origin << "\tr.direction: " << r.direction << endl
                   << "rad: " << rad << "\tcenter: " << center << endl);
    float a, b, c, disc;
    float t1, t2;
    a = Sqr(r.direction.x) + Sqr(r.direction.y);
    b = 2.0 * ((r.origin.x - center.x) * r.direction.x + (r.origin.y - center.y) * r.direction.y);
    c = Sqr(r.origin.x - center.x) + Sqr(r.origin.y - center.y) - Sqr(rad);
    DebugGeom(printf(" RCI: a: %f\tb: %f\t c: %f\n", a, b, c));

    disc = Sqr(b) - 4 * a * c;
    if (disc < 0)
    {
        DebugGeom(printf(" RCI disc < 0: %f\n", disc));
        return 0;
    }

    disc = sqrt(disc);
    t1 = (-b + disc) / (2.0 * a);
    t2 = (-b - disc) / (2.0 * a);
    DebugGeom(printf(" RCI: t1: %f\tt2: %f\n", t1, t2));

    if (t1 > t2)
    {
        DebugGeom(printf(" RCI: reversing t1, t2\n"));
        float temp = t1;
        t1 = t2;
        t2 = temp;
    }

    if (t1 > 0.0)
    {
        if (t2 > 0.0)
        {
            *psol1 = r.origin + r.direction * t1;
            *psol2 = r.origin + r.direction * t2;
            DebugGeom(printf(" RCI:two sols\n"));
            return 2;
        }
        else
        {
            my_error("RayCircleIntersect: weird roots");
            return 0;
        }
    }
    else if (t2 > 0.0)
    {
        *psol1 = r.origin + r.direction * t2;
        DebugGeom(printf(" RCI:t2 only sol\n"));
        return 1;
    }
    else
        return 0;

    return 0;
}
#endif

/********************************************************************************/

int QuadraticFormula(float a, float b, float c, float *psol1, float *psol2)
{
    float d = Sqr(b) - 4 * a * c;
    if (fabs(d) < FLOAT_EPS)
    {
        *psol1 = -b / (2 * a);
        return 1;
    }
    else if (d < 0)
    {
        return 0;
    }
    else
    {
        d = sqrt(d);
        *psol1 = (-b + d) / (2 * a);
        *psol2 = (-b - d) / (2 * a);
        return 2;
    }
}

/********************************************************************************/

/* some test code
   Vector sol1, sol2;
  int num;
  num = LineCircleIntersect(LineFromTwoPoints(Vector(-10, 2), Vector(10, 2)),
                          1, Vector(1,2), &sol1, &sol2);
  cout << "num: " << num << "\tsol1: " << sol1 << "\tsol2: " << sol2 << endl;
  num = LineCircleIntersect(LineFromTwoPoints(Vector(2, 10), Vector(2, -10)),
                          1, Vector(2,1), &sol1, &sol2);
  cout << "num: " << num << "\tsol1: " << sol1 << "\tsol2: " << sol2 << endl;
  num = LineCircleIntersect(LineFromTwoPoints(Vector(-10, -10), Vector(10, 10)),
                          sqrt(2), Vector(1,1), &sol1, &sol2);
  cout << "num: " << num << "\tsol1: " << sol1 << "\tsol2: " << sol2 << endl;
  exit(1);
  */
int LineCircleIntersect(Line l, float rad, Vector center,
                        Vector *psol1, Vector *psol2)
{
    *psol1 = *psol2 = Vector(0, 0);
    if (fabs(l.A) > FLOAT_EPS)
    {
        float a, b, c;
        a = 1 + Sqr(l.B / l.A);
        b = 2 * (-center.y + (l.C / l.A + center.x) * (l.B / l.A));
        c = -Sqr(rad) + Sqr(l.C / l.A + center.x) + Sqr(center.y);
        int numSol = QuadraticFormula(a, b, c, &(psol1->y), &(psol2->y));
        psol1->x = -((l.B * psol1->y + l.C) / l.A);
        psol2->x = -((l.B * psol2->y + l.C) / l.A);
        return numSol;
    }
    else
    {
        int numSol = QuadraticFormula(1, -2 * center.x,
                                      Sqr(center.x) + Sqr(l.C / l.B + center.y) - Sqr(rad),
                                      &(psol1->x), &(psol2->x));
        psol1->y = psol2->y = -l.C / l.B;
        return numSol;
    }
}

/********************************************************************************/

/* loop over rectangle edges. See if we're outside of the edge.
   If so, try to intersect with edge to move it inside */
Vector AdjustPtToRectOnLine(Vector pt, Rectangle r, Line l)
{
    DebugGeom(cout << "Adjust: " << pt << "\t" << r << "\t" << l << endl);
    Vector c = r.Center();
    for (int i = 0; i < 4; i++)
    {
        Line edge = r.GetEdge(i);
        if (edge.HalfPlaneTest(pt) != edge.HalfPlaneTest(c) && !edge.SameSlope(l))
        {
            // pt is outside this edge
            DebugGeom(printf("AdjustPtToRectOnLine: HalfPlaneTest failed for %d\n", i));
            Vector newPt = edge.intersection(l);
            if (edge.InBetween(newPt, r.GetPoint(i), r.GetPoint(i + 1)))
            {
                DebugGeom(printf("AdjustPtToRectOnLine: Intersection okay for %d\n", i));
                return newPt;
            }
        }
    }
    return pt;
}

/********************************************************************************/

Bool InBetween(Vector pt, Vector end1, Vector end2)
{
    Line l = LineFromTwoPoints(end1, end2);
    return l.InBetween(pt, end1, end2);
}

/********************************************************************************/

Vector PointInBetween(Vector pt1, Vector pt2, float pt1dist)
{
    if (pt1dist < 0)
        my_error("no neg dists");
    Vector targ = pt2 - pt1;
    targ *= pt1dist / targ.mod();
    return (pt1 + targ);
}

/********************************************************************************/

AngleDeg AngleBisect(AngleDeg a1, AngleDeg a2)
{
    if (fabs(a1 - a2) > 180)
        return GetNormalizeAngleDeg((Min(a1, a2) + 360 + Max(a1, a2)) / 2);

    return GetNormalizeAngleDeg((a1 + a2) / 2);
}

/********************************************************************************/

Vector GetClosestPtInBetween(Vector pt, Vector end1, Vector end2)
{
    Line l;
    l.LineFromTwoPoints(end1, end2);
    return l.GetClosestPtInBetween(pt, end1, end2);
}

/********************************************************************************/

Bool IsPointInCone(Vector pt, float wid_dist_ratio, Vector end, Vector vert)
{
    Line l = LineFromTwoPoints(vert, end);
    Vector proj_pt = l.ProjectPoint(pt);
    return ((proj_pt.dist2(pt) < proj_pt.dist2(vert) * wid_dist_ratio * wid_dist_ratio &&
             l.InBetween(proj_pt, vert, end)))
               ? TRUE
               : FALSE;
}

/* -*- Mode: C++ -*- */

/* MemOption.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

/* setting defaults to match version 4.06 server.conf */
OptionInfo::OptionInfo()
{

    VP_test_l = FALSE;
    VP_test_r = FALSE;
    VP_test = FALSE;
    VP_train_DT = FALSE;
    VP_use_DT = FALSE;

    IP_my_score = 0;
    IP_their_score = 0;
    IP_reconnect = 0;

    sprintf(MyTeamName, "XMUnited");

    /* no option flags for these */
    SP_pitch_length = 105.0;
    SP_pitch_width = 68.0;
    SP_pitch_margin = 5.0;
    SP_penalty_area_length = 16.5;
    SP_penalty_area_width = 40.32;
    SP_goal_area_length = 5.5;
    SP_goal_area_width = 18.32;
    SP_penalty_spot_dist = 11.0;
    SP_corner_arc_r = 1.0;
    SP_free_kick_buffer = 9.15;
    SP_after_goal_wait = 50;
    SP_feel_distance = 3.0;
    SP_num_lines = 4;
    SP_num_markers = 55;
    SP_unum_far_length = 20.0;
    SP_unum_too_far_length = 40.0;
    SP_team_far_length = 40.0;
    SP_team_too_far_length = 60.0;

    SP_version = 5.18;
    SP_team_size = 11;
    SP_half = 1;
    sprintf(SP_host, "127.0.0.1");
    SP_goal_width = 14.02;
    SP_player_size = 0.8;
    SP_player_decay = 0.4;
    SP_player_rand = 0.1;
    SP_player_weight = 60.0;
    SP_player_speed_max = 32.0;
    SP_stamina_max = 2000.0;
    SP_stamina_inc = 20.0;
    SP_recover_dec_thr = 0.3;
    SP_recover_dec = 0.05;
    SP_recover_min = 0.1;
    SP_effort_dec_thr = 0.4;
    SP_effort_min = 0.1;
    SP_effort_dec = 0.05;
    SP_effort_inc_thr = 0.9;
    SP_effort_inc = 0.05;
    SP_ball_size = 0.085;
    SP_ball_decay = 0.96;
    SP_ball_rand = 0.05;
    SP_ball_weight = 0.2;
    ;
    SP_ball_speed_max = 32.0;
    SP_dash_power_rate = 0.01;
    SP_kick_power_rate = 0.016;
    SP_kickable_margin = 1.0;
    SP_kickable_area = SP_kickable_margin + SP_ball_size + SP_player_size;
    SP_catch_prob = 1.0;
    SP_catch_area_l = 2.0;
    SP_catch_area_w = 1.0;
    SP_catch_ban_cycle = 5;
    SP_max_power = 100;
    SP_min_power = -30;
    SP_max_moment = 180;
    SP_min_moment = -180;
    SP_min_neck_angle = -90.0;
    SP_max_neck_angle = 90.0;
    SP_min_neck_moment = -180.0;
    SP_max_neck_moment = 180.0;
    SP_visible_angle = 90.0;
    SP_audio_cut_dist = 50.0;
    SP_dist_qstep = 0.1;
    SP_land_qstep = 0.01;
    SP_ckmargin = 1.0;
    SP_wind_dir = 0.0;
    SP_wind_force = 0.0;
    SP_wind_rand = 0.0;
    SP_wind_none = FALSE;
    SP_wind_random = FALSE;
    SP_half_time = 300;
    SP_port = 6000;
    SP_coach_port = 6001;
    SP_olcoach_port = 6002;
    SP_simulator_step = 100;
    SP_send_step = 150;
    SP_recv_step = 10;
    SP_say_msg_size = 512;
    SP_hear_max = 2;
    SP_hear_inc = 1;
    SP_hear_decay = 2;
    SP_coach_mode = FALSE;
    SP_coach_w_referee_mode = FALSE;
    SP_say_coach_cnt_max = 128;
    SP_say_coach_msg_size = 128;
    SP_send_vi_step = 100;
    SP_look_step = 100;
    SP_use_offside = FALSE;
    SP_forbid_kickoff_offside = TRUE;
    SP_verbose = TRUE;
    SP_offside_area = 9.15;
    SP_inertia_moment = 5.0;
    SP_sense_body_step = 100;
    SP_offside_kick_margin = 9.15;
    SP_record_messages = FALSE;

    CP_goalie = FALSE;
    CP_save_log = FALSE;
    CP_save_freq = 10;
    CP_save_sound_log = FALSE;
    CP_save_sound_freq = 10;
    CP_save_action_log_level = 0; /* 0 means save nothing */
    CP_save_action_freq = 40;
    CP_send_ban_recv_step_factor = 3.0;
    CP_interrupts_per_cycle = 5;
    CP_interrupts_left_to_act = 5;
    CP_max_conf = 1.0;
    CP_min_valid_conf = 0.5;
    CP_conf_decay = 0.98;
    CP_player_conf_decay = 0.99;
    CP_ball_conf_decay = 0.95;
    CP_max_player_move_factor = 7;
    CP_max_say_interval = 100;
    CP_ball_moving_threshold = .2; /* experimentally checked -- ball still, player fast => .15 ball speed */
    CP_dodge_angle_buffer = 25;
    CP_dodge_distance_buffer = 3.5;
    CP_dodge_power = 100;
    CP_dodge_angle = 60;
    sprintf(CP_tree_stem, "pass");
    CP_DT_evaluate_interval = 10;
    CP_say_tired_interval = 20;
    CP_tired_buffer = 10;
    CP_set_plays = TRUE;
    CP_Setplay_Delay = 5;
    CP_Setplay_Say_Delay = SP_hear_decay * 5;
    CP_Setplay_Max_Delay = 100;
    CP_Setplay_Time_Limit = 150;
    CP_kickable_buffer = .2;
    CP_mark_persist_time = 6000;
    CP_track_min_distance = 0;
    CP_track_max_distance = 15;
    CP_pull_offsides = FALSE;
    CP_pull_offsides_when_winning = TRUE;
    CP_spar = TRUE;
    CP_mark = TRUE;
    CP_communicate = FALSE;
    CP_change_view_for_ball_cycles = 2;
    CP_defer_kick_to_teammate_buffer = .05;
    CP_scan_overlap_angle = 2;

    CP_pull_offside_threshold = 5;
    CP_pull_offside_buffer = 3;

    CP_ball_forget_angle_buf = 3;
    CP_player_forget_angle_buf = 5;
    CP_ball_forget_dist_buf = 1;
    CP_player_forget_dist_buf = 1;

    CP_beat_offsides_buffer = 15;
    CP_beat_offsides_threshold = 30;
    CP_beat_offsides_max_x = 25;
    CP_congestion_epsilon = .01;
    CP_back_pass_opponent_buffer = 10;
    CP_back_pass_offside_buffer = 10;
    CP_min_less_congested_pass_dist = 7;

    /* pat added these */
    CP_use_new_position_based_vel = TRUE;
    CP_stop_on_error = FALSE;

    CP_opt_ctrl_dist = .7; //(SP_player_size + .8 *SP_kickable_margin);
    CP_KickTo_err = 3;
    CP_closest_margin = .585; //(SP_player_size + 2.0*SP_ball_size);
    CP_dokick_factor = .4;
    CP_max_turn_kick_pow = 30;

    CP_max_ignore_vel = .0;
    CP_kick_time_space = 1;
    CP_max_est_err = .3;
    CP_max_go_to_point_angle_err = 5;
    CP_holdball_kickable_buffer = .1;
    CP_stop_ball_power = 30;
    CP_possessor_intercept_space = 4;
    CP_can_keep_ball_cycle_buffer = 1;

    /* no longer used
    CP_hard_kick_margin     = .97;//(SP_player_size + 2.0*SP_ball_size);
    CP_hard_kick_factor     = .25;
    CP_hard_kick_end_turn_dist = 1.1;//(SP_player_size + .3 *SP_kickable_margin);
    */
    CP_hard_kick_dist_buffer = .1;
    CP_max_hard_kick_angle_err = 5;

    CP_hardest_kick_player_ang = 90;  // angle relative to direction of ball
    CP_hardest_kick_ball_dist = .93; // kickable_area * .6
    CP_hardest_kick_ball_ang = 180;    // this is realtive to the direction of travel
    CP_max_dash_help_kick_angle = 60;

    CP_max_int_lookahead = 50;
    CP_intercept_step = 5;
    CP_my_intercept_step = 1;
    CP_intercept_aim_ahead = 1;
    CP_no_turn_max_cyc_diff = 2;
    CP_no_turn_max_dist_diff = 1.0;
    CP_turnball_opp_worry_dist = 5;
    CP_collision_buffer = .2;
    CP_behind_angle = 80;
    CP_time_for_full_rotation = 12; /* average guestimate */
    CP_ball_vel_invalidation_factor = 2.0;

    CP_dribble_dash_pow = 100;
    CP_dribble_ball_dist = .75;
    CP_dribble_angle_norm = 60;
    CP_dribble_exp_angle_buffer = 10;

    CP_dribble_ignore_opp_dist = 15;
    CP_dribble_worry_opp_dist = 4;
    CP_dribble_dodge_max_dist = 6;
    CP_dribble_dodge_angle_err = 15;
    CP_dribble_angle_ignore_buffer = 5;
    CP_dribble_dodge_close_dist = 2;
    CP_dribble_scan_field = TRUE;
    CP_can_dribble_cone_ratio = .5;

    CP_dribble_towards_length = 10;
    CP_dribble_sideline_buffer = 1.5;
    CP_dribble_circle_inner_rad = 2.5;
    CP_dribble_circle_outer_rad = 4;
    CP_dribble_circle_ang = 90;

    CP_move_imp_1v1_initial = 0.0;
    CP_move_imp_1v1_inc = .2;
    CP_move_imp_1v1_threshold = 1.0;
    CP_at_point_buffer = 1;
    CP_overrun_dist = 3;
    CP_def_block_dist = 2;
    CP_def_block_dist_ratio = .5;
    CP_overrun_buffer = 2.5;
    CP_cycles_to_kick = 3;
    CP_breakaway_buffer = 2;
    CP_our_breakaway_kickable_buffer = 1.5;
    CP_their_breakaway_front_kickable_buffer = 5.0;
    CP_their_breakaway_back_kickable_buffer = 2.0;

    CP_breakaway_approach_x = 35;
    CP_breakaway_approach_y = 8;
    CP_breakaway_targ_valid_time = 3;
    CP_breakaway_min_goalie_steal_time = 6;
    CP_breakaway_kick_run_min_cycles = 7;
    CP_breakaway_kick_run_max_cycles = 2;
    CP_our_breakaway_min_cone_dist_wid = 18;
    CP_their_breakaway_min_cone_dist_wid = 12;
    CP_breakaway_middle_buffer = 3;
    CP_breakaway_kick_run_worry_dist = 10;
    CP_goalie_breakaway_kickable_buffer = 1.5;
    CP_breakaway_mode = 0;

    CP_static_kick_dist_err = .16; // old: .14
    CP_static_kick_ang_err = 4;  // old: 5
    // no longer used
    // CP_static_kick_dist = .985;
    // CP_static_kick_ang = 47;  /* caculated value! */
    // CP_static_kick_ang = 42;  /* caculated value!- extar buffer */
    //   CP_static_kick_overrun_dist = 4;

    CP_goalie_baseline_buffer = 1;
    CP_goalie_scan_angle_err = 10;
    CP_goalie_at_point_buffer = .1;
    CP_goalie_vis_angle_err = 5;
    CP_goalie_max_shot_distance = 40;
    CP_goalie_min_pos_dist = 15;
    CP_goalie_max_pos_dist = SP_pitch_length * .75;
    CP_goalie_max_forward_percent = .75;
    CP_goalie_ball_ang_for_corner = 90;
    CP_goalie_max_come_out_dist = 10;
    CP_goalie_ball_dist_for_corner = SP_penalty_area_length;
    CP_goalie_ball_dist_for_center = SP_pitch_length / 2;
    CP_goalie_free_kick_dist = 3;
    CP_goalie_go_to_ball_cone_ratio = .25;
    CP_goalie_warn_space = 10;
    // CP_goalie_comes_out = TRUE;
    CP_goalie_comes_out = FALSE;
    CP_goalie_catch_wait_time = 2;
    CP_goalie_opponent_dist_to_block = 24;
    CP_goalie_position_weight_dist = 10;
    CP_goalie_narrow_sideline_cyc = 3;
    CP_goalie_no_buffer_dist = 10;

    CP_clear_ball_ang_step = 5.0;
    CP_clear_ball_cone_ratio = .5;
    CP_clear_ball_max_dist = 30;
    CP_clear_offensive_min_horiz_dist = 20;
    CP_clear_offensive_min_angle = 60;

    CP_should_cross_corner_dist = 10; // pitch_width /2 - penalty_area_w / 2
    CP_should_cross_baseline_buffer = 6;
    CP_should_move_to_cross_corner_dist = 20;
    CP_cross_pt_x = 36; // pitch_length / 2 - penalty_area_l
    CP_cross_pt_y = 9;  // goalie_area_w / 2
    CP_cross_target_vel = .5;

    CP_dont_dribble_to_middle_min_x = 20;

    /* not used anymore
       CP_hardest_kick_shot_distance = 13;
       CP_moderate_kick_shot_distance = 9;
    */
    CP_good_shot_distance = 17;
    CP_shot_distance = 25;
    CP_cycles_to_kick_buffer = 3;
    CP_better_shot_cyc_diff = 5;
    // CP_breakaway_shot_distance = 16;
    CP_shot_speed = 1.7; // our average shot speed
    CP_shot_goalie_react_buffer = 6;
    CP_good_shot_goalie_react_buffer = 3;

    sprintf(FP_initial_formation, "433");
    sprintf(FP_formation_when_tied, "433");
    sprintf(FP_formation_when_losing, "334");
    sprintf(FP_formation_when_losing_lots, "334");
    sprintf(FP_formation_when_winning, "532");
    sprintf(FP_initial_hc_method, "Shift");
    sprintf(FP_initial_mc_method, "Obey");
    FP_initial_player_1_pos = 1;
    FP_initial_player_2_pos = 2;
    FP_initial_player_3_pos = 3;
    FP_initial_player_4_pos = 4;
    FP_initial_player_5_pos = 5;
    FP_initial_player_6_pos = 6;
    FP_initial_player_7_pos = 7;
    FP_initial_player_8_pos = 8;
    FP_initial_player_9_pos = 9;
    FP_initial_player_10_pos = 10;
    FP_initial_player_11_pos = 11;
    FP_goalie_number = 11;
}

void OptionInfo::GetOptions(int argc, char **argv)
{
    option_t opt[] = {
        {"test_l", (void *)&VP_test_l, V_ONOFF},
        {"test_r", (void *)&VP_test_r, V_ONOFF},
        {"test", (void *)&VP_test, V_ONOFF},
        {"train_DT", (void *)&VP_train_DT, V_ONOFF},
        {"use_DT", (void *)&VP_use_DT, V_ONOFF},

        {"my_score", (void *)&IP_my_score, V_INT},
        {"their_score", (void *)&IP_their_score, V_INT},
        {"reconnect", (void *)&IP_reconnect, V_INT},

        {"team_name", (void *)&MyTeamName, V_STRING},
        {"goalie", (void *)&CP_goalie, V_ONOFF},
        {"save_log", (void *)&CP_save_log, V_ONOFF},
        {"save_freq", (void *)&CP_save_freq, V_INT},
        {"save_sound_log", (void *)&CP_save_sound_log, V_ONOFF},
        {"save_sound_freq", (void *)&CP_save_sound_freq, V_INT},
        {"save_action_log_level", (void *)&CP_save_action_log_level, V_INT},
        {"save_action_freq", (void *)&CP_save_action_freq, V_INT},
        {"send_ban_recv_step_factor", (void *)&CP_send_ban_recv_step_factor, V_FLOAT},
        {"interrupts_per_cycle", (void *)&CP_interrupts_per_cycle, V_INT},
        {"interrupts_left_to_act", (void *)&CP_interrupts_left_to_act, V_INT},
        {"max_conf", (void *)&CP_max_conf, V_FLOAT},
        {"min_conf", (void *)&CP_min_valid_conf, V_FLOAT},
        {"conf_decay", (void *)&CP_conf_decay, V_FLOAT},
        {"player_conf_decay", (void *)&CP_player_conf_decay, V_FLOAT},
        {"ball_conf_decay", (void *)&CP_ball_conf_decay, V_FLOAT},
        {"max_player_move_factor", (void *)&CP_max_player_move_factor, V_FLOAT},
        {"max_say_interval", (void *)&CP_max_say_interval, V_INT},
        {"ball_moving_threshold", (void *)&CP_ball_moving_threshold, V_FLOAT},
        {"dodge_distance_buffer", (void *)&CP_dodge_distance_buffer, V_FLOAT},
        {"dodge_angle_buffer", (void *)&CP_dodge_angle_buffer, V_FLOAT},
        {"dodge_power", (void *)&CP_dodge_power, V_FLOAT},
        {"dodge_angle", (void *)&CP_dodge_angle, V_FLOAT},
        {"tree_stem", (void *)&CP_tree_stem, V_STRING},
        {"DT_evaluate_interval", (void *)&CP_DT_evaluate_interval, V_INT},
        {"say_tired_interval", (void *)&CP_say_tired_interval, V_INT},
        {"tired_buffer", (void *)&CP_tired_buffer, V_FLOAT},
        {"set_plays", (void *)&CP_set_plays, V_ONOFF},
        {"set_play_delay", (void *)&CP_Setplay_Delay, V_INT},
        {"set_play_say_delay", (void *)&CP_Setplay_Say_Delay, V_INT},
        {"set_play_time_limit", (void *)&CP_Setplay_Time_Limit, V_INT},
        {"kickable_buffer", (void *)&CP_kickable_buffer, V_FLOAT},
        {"mark_persist_time", (void *)&CP_mark_persist_time, V_INT},
        {"track_max_distance", (void *)&CP_track_max_distance, V_FLOAT},
        {"track_min_distance", (void *)&CP_track_min_distance, V_FLOAT},
        {"pull_offsides", (void *)&CP_pull_offsides, V_ONOFF},
        {"pull_offsides_when_winning", (void *)&CP_pull_offsides_when_winning, V_ONOFF},
        {"spar", (void *)&CP_spar, V_ONOFF},
        {"mark", (void *)&CP_mark, V_ONOFF},
        {"communicate", (void *)&CP_communicate, V_ONOFF},
        {"change_view_for_ball_cycles", (void *)&CP_change_view_for_ball_cycles, V_INT},
        {"defer_kick_to_teammate_buffer", (void *)&CP_defer_kick_to_teammate_buffer, V_FLOAT},
        {"scan_overlap_angle", (void *)&CP_scan_overlap_angle, V_FLOAT},

        {"pull_offside_threshold", (void *)&CP_pull_offside_threshold, V_FLOAT},
        {"pull_offside_buffer", (void *)&CP_pull_offside_buffer, V_FLOAT},

        {"ball_forget_angle_buf", (void *)&CP_ball_forget_angle_buf, V_FLOAT},
        {"player_forget_angle_buf", (void *)&CP_player_forget_angle_buf, V_FLOAT},
        {"ball_forget_dist_buf", (void *)&CP_ball_forget_dist_buf, V_FLOAT},
        {"player_forget_dist_buf", (void *)&CP_player_forget_dist_buf, V_FLOAT},

        {"beat_offsides_buffer", (void *)&CP_beat_offsides_buffer, V_FLOAT},
        {"beat_offsides_threshold", (void *)&CP_beat_offsides_threshold, V_FLOAT},
        {"beat_offsides_max_x", (void *)&CP_beat_offsides_max_x, V_FLOAT},
        {"congestion_epsilon", (void *)&CP_congestion_epsilon, V_FLOAT},
        {"back_pass_opponent_buffer", (void *)&CP_back_pass_opponent_buffer, V_FLOAT},
        {"back_pass_offside_buffer", (void *)&CP_back_pass_offside_buffer, V_FLOAT},
        {"min_less_congested_pass_dist", (void *)&CP_min_less_congested_pass_dist, V_FLOAT},

        {"use_new_position_based_vel", (void *)&CP_use_new_position_based_vel, V_ONOFF},
        {"stop_on_error", (void *)&CP_stop_on_error, V_ONOFF},

        {"opt_ctrl_dist", (void *)&CP_opt_ctrl_dist, V_FLOAT},
        {"KickTo_err", (void *)&CP_KickTo_err, V_FLOAT},
        {"closest_margin", (void *)&CP_closest_margin, V_FLOAT},
        {"dokick_factor", (void *)&CP_dokick_factor, V_FLOAT},
        {"max_turn_kick_pow", (void *)&CP_max_turn_kick_pow, V_FLOAT},
        {"kick_time_space", (void *)&CP_kick_time_space, V_INT},
        {"max_ignore_vel", (void *)&CP_max_ignore_vel, V_FLOAT},
        {"max_est_err", (void *)&CP_max_est_err, V_FLOAT},
        {"holdball_kickable_buffer", (void *)&CP_holdball_kickable_buffer, V_FLOAT},
        {"stop_ball_power", (void *)&CP_stop_ball_power, V_INT},
        {"possessor_intercept_space", (void *)&CP_possessor_intercept_space, V_INT},
        {"can_keep_ball_cycle_buffer", (void *)&CP_can_keep_ball_cycle_buffer, V_INT},

        /* no longer used
        {"hard_kick_margin",            (void *)&CP_hard_kick_margin,   V_FLOAT},
        {"hard_kick_end_turn_dist",     (void *)&CP_hard_kick_end_turn_dist,  V_FLOAT},
        {"hard_kick_factor",            (void *)&CP_hard_kick_factor,   V_FLOAT},
        */
        {"max_hard_kick_angle_err", (void *)&CP_max_hard_kick_angle_err, V_INT},
        {"hard_kick_dist_buffer", (void *)&CP_hard_kick_dist_buffer, V_FLOAT},
        {"hardest_kick_ball_ang", (void *)&CP_hardest_kick_ball_ang, V_INT},
        {"hardest_kick_ball_dist", (void *)&CP_hardest_kick_ball_dist, V_FLOAT},
        {"hardest_kick_player_ang", (void *)&CP_hardest_kick_player_ang, V_INT},
        {"max_dash_help_kick_angle", (void *)&CP_max_dash_help_kick_angle, V_FLOAT},

        {"max_go_to_point_angle_err", (void *)&CP_max_go_to_point_angle_err, V_INT},
        {"max_int_lookahead", (void *)&CP_max_int_lookahead, V_INT},
        {"intercept_close_dist", (void *)&CP_intercept_close_dist, V_FLOAT},
        {"intercept_step", (void *)&CP_intercept_step, V_INT},
        {"my_intercept_step", (void *)&CP_my_intercept_step, V_INT},
        {"intercept_aim_ahead", (void *)&CP_intercept_aim_ahead, V_INT},
        {"no_turn_max_cyc_diff", (void *)&CP_no_turn_max_cyc_diff, V_INT},
        {"no_turn_max_dist_diff", (void *)&CP_no_turn_max_dist_diff, V_FLOAT},

        {"turnball_opp_worry_dist", (void *)&CP_turnball_opp_worry_dist, V_FLOAT},
        {"collision_buffer", (void *)&CP_collision_buffer, V_FLOAT},
        {"behind_angle", (void *)&CP_behind_angle, V_FLOAT},
        {"ball_vel_invalidation_factor", (void *)&CP_ball_vel_invalidation_factor, V_FLOAT},
        {"time_for_full_rotation", (void *)&CP_time_for_full_rotation, V_INT},

        {"dribble_dash_pow", (void *)&CP_dribble_dash_pow, V_INT},
        {"dribble_ball_dist", (void *)&CP_dribble_ball_dist, V_FLOAT},
        {"dribble_ignore_opp_dist", (void *)&CP_dribble_ignore_opp_dist, V_FLOAT},
        {"dribble_worry_opp_dist", (void *)&CP_dribble_worry_opp_dist, V_FLOAT},
        {"dribble_angle_norm", (void *)&CP_dribble_angle_norm, V_FLOAT},
        {"dribble_dodge_max_dist", (void *)&CP_dribble_dodge_max_dist, V_FLOAT},
        {"dribble_dodge_angle_err", (void *)&CP_dribble_dodge_angle_err, V_FLOAT},
        {"dribble_exp_angle_buffer", (void *)&CP_dribble_exp_angle_buffer, V_FLOAT},
        {"dribble_angle_ignore_buffer", (void *)&CP_dribble_angle_ignore_buffer, V_FLOAT},
        {"dribble_dodge_close_dist", (void *)&CP_dribble_dodge_close_dist, V_FLOAT},
        {"dribble_scan_field", (void *)&CP_dribble_scan_field, V_ONOFF},

        {"can_dribble_cone_ratio", (void *)&CP_can_dribble_cone_ratio, V_FLOAT},
        {"dribble_towards_length", (void *)&CP_dribble_towards_length, V_FLOAT},
        {"dribble_sideline_buffer", (void *)&CP_dribble_sideline_buffer, V_FLOAT},
        {"dribble_circle_inner_rad", (void *)&CP_dribble_circle_inner_rad, V_FLOAT},
        {"dribble_circle_outer_rad", (void *)&CP_dribble_circle_outer_rad, V_FLOAT},
        {"dribble_circle_ang", (void *)&CP_dribble_circle_ang, V_FLOAT},

        {"move_imp_1v1_initial", (void *)&CP_move_imp_1v1_initial, V_FLOAT},
        {"move_imp_1v1_inc", (void *)&CP_move_imp_1v1_inc, V_FLOAT},
        {"move_imp_1v1_threshold", (void *)&CP_move_imp_1v1_threshold, V_FLOAT},
        {"at_point_buffer", (void *)&CP_at_point_buffer, V_FLOAT},
        {"overrun_dist", (void *)&CP_overrun_dist, V_FLOAT},
        {"def_block_dist", (void *)&CP_def_block_dist, V_FLOAT},
        {"def_block_dist_ratio", (void *)&CP_def_block_dist_ratio, V_FLOAT},
        {"overrun_buffer", (void *)&CP_overrun_buffer, V_FLOAT},
        {"cycles_to_kick", (void *)&CP_cycles_to_kick, V_FLOAT},
        {"breakaway_buffer", (void *)&CP_breakaway_buffer, V_FLOAT},
        {"our_breakaway_kickable_buffer", (void *)&CP_our_breakaway_kickable_buffer, V_FLOAT},
        {"their_breakaway_front_kickable_buffer", (void *)&CP_their_breakaway_front_kickable_buffer, V_FLOAT},
        {"their_breakaway_back_kickable_buffer", (void *)&CP_their_breakaway_back_kickable_buffer, V_FLOAT},
        {"goalie_breakaway_kickable_buffer", (void *)&CP_goalie_breakaway_kickable_buffer, V_FLOAT},

        {"breakaway_approach_x", (void *)&CP_breakaway_approach_x, V_FLOAT},
        {"breakaway_approach_y", (void *)&CP_breakaway_approach_y, V_FLOAT},
        {"breakaway_targ_valid_time", (void *)&CP_breakaway_targ_valid_time, V_INT},
        {"breakaway_min_goalie_steal_time", (void *)&CP_breakaway_min_goalie_steal_time, V_INT},
        {"breakaway_kick_run_min_cycles", (void *)&CP_breakaway_kick_run_min_cycles, V_INT},
        {"breakaway_kick_run_max_cycles", (void *)&CP_breakaway_kick_run_max_cycles, V_INT},
        {"their_breakaway_min_cone_dist_wid", (void *)&CP_their_breakaway_min_cone_dist_wid, V_FLOAT},
        {"our_breakaway_min_cone_dist_wid", (void *)&CP_our_breakaway_min_cone_dist_wid, V_FLOAT},
        {"breakaway_middle_buffer", (void *)&CP_breakaway_middle_buffer, V_FLOAT},
        {"breakaway_kick_run_worry_dist", (void *)&CP_breakaway_kick_run_worry_dist, V_FLOAT},
        {"breakaway_mode", (void *)&CP_breakaway_mode, V_INT},

        {"static_kick_dist_err", (void *)&CP_static_kick_dist_err, V_FLOAT},
        {"static_kick_ang_err", (void *)&CP_static_kick_ang_err, V_FLOAT},
        /* no longer used
        {"static_kick_dist",            (void *)&CP_static_kick_dist, V_FLOAT},
        {"static_kick_ang",             (void *)&CP_static_kick_ang, V_FLOAT},
        {"static_kick_overrun_dist",    (void *)&CP_static_kick_overrun_dist, V_FLOAT},
        */

        {"goalie_baseline_buffer", (void *)&CP_goalie_baseline_buffer, V_FLOAT},
        {"goalie_scan_angle_err", (void *)&CP_goalie_scan_angle_err, V_FLOAT},
        {"goalie_at_point_buffer", (void *)&CP_goalie_at_point_buffer, V_FLOAT},
        {"goalie_vis_angle_err", (void *)&CP_goalie_vis_angle_err, V_FLOAT},
        {"goalie_max_shot_distance", (void *)&CP_goalie_max_shot_distance, V_FLOAT},
        {"goalie_min_pos_dist", (void *)&CP_goalie_min_pos_dist, V_FLOAT},
        {"goalie_max_pos_dist", (void *)&CP_goalie_max_pos_dist, V_FLOAT},
        {"goalie_max_forward_percent", (void *)&CP_goalie_max_forward_percent, V_FLOAT},
        {"goalie_ball_ang_for_corner", (void *)&CP_goalie_ball_ang_for_corner, V_FLOAT},
        {"goalie_max_come_out_dist", (void *)&CP_goalie_max_come_out_dist, V_FLOAT},
        {"goalie_ball_dist_for_corner", (void *)&CP_goalie_ball_dist_for_corner, V_FLOAT},
        {"goalie_ball_dist_for_center", (void *)&CP_goalie_ball_dist_for_center, V_FLOAT},
        {"goalie_free_kick_dist", (void *)&CP_goalie_free_kick_dist, V_FLOAT},
        {"goalie_go_to_ball_cone_ratio", (void *)&CP_goalie_go_to_ball_cone_ratio, V_FLOAT},
        {"goalie_warn_space", (void *)&CP_goalie_warn_space, V_INT},
        {"goalie_comes_out", (void *)&CP_goalie_comes_out, V_ONOFF},
        {"goalie_catch_wait_time", (void *)&CP_goalie_catch_wait_time, V_INT},
        {"goalie_opponent_dist_to_block", (void *)&CP_goalie_opponent_dist_to_block, V_FLOAT},
        {"goalie_position_weight_dist", (void *)&CP_goalie_position_weight_dist, V_FLOAT},
        {"goalie_narrow_sideline_cyc", (void *)&CP_goalie_narrow_sideline_cyc, V_INT},
        {"goalie_no_buffer_dist", (void *)&CP_goalie_no_buffer_dist, V_FLOAT},

        {"clear_ball_ang_step", (void *)&CP_clear_ball_ang_step, V_FLOAT},
        {"clear_ball_cone_ratio", (void *)&CP_clear_ball_cone_ratio, V_FLOAT},
        {"clear_ball_max_dist", (void *)&CP_clear_ball_max_dist, V_FLOAT},
        {"clear_offensive_min_horiz_dist", (void *)&CP_clear_offensive_min_horiz_dist, V_FLOAT},
        {"clear_offensive_min_angle", (void *)&CP_clear_offensive_min_angle, V_FLOAT},

        {"should_cross_corner_dist", (void *)&CP_should_cross_corner_dist, V_FLOAT},
        {"should_cross_baseline_buffer", (void *)&CP_should_cross_baseline_buffer, V_FLOAT},
        {"should_move_to_cross_corner_dist", (void *)&CP_should_move_to_cross_corner_dist, V_FLOAT},
        {"cross_pt_x", (void *)&CP_cross_pt_x, V_FLOAT},
        {"cross_pt_y", (void *)&CP_cross_pt_y, V_FLOAT},
        {"cross_target_vel", (void *)&CP_cross_target_vel, V_FLOAT},

        {"dont_dribble_to_middle_min_x", (void *)&CP_dont_dribble_to_middle_min_x, V_FLOAT},

        /* not used anymore
          {"hardest_kick_shot_distance",  (void *)&CP_hardest_kick_shot_distance, V_FLOAT},
          {"moderate_kick_shot_distance", (void *)&CP_moderate_kick_shot_distance, V_FLOAT},
          */
        {"good_shot_distance", (void *)&CP_good_shot_distance, V_FLOAT},
        {"shot_distance", (void *)&CP_shot_distance, V_FLOAT},
        {"cycles_to_kick_buffer", (void *)&CP_cycles_to_kick_buffer, V_INT},
        {"shot_speed", (void *)&CP_shot_speed, V_FLOAT},
        {"shot_goalie_react_buffer", (void *)&CP_shot_goalie_react_buffer, V_INT},
        {"good_shot_goalie_react_buffer", (void *)&CP_good_shot_goalie_react_buffer, V_INT},
        {"better_shot_cyc_diff", (void *)&CP_better_shot_cyc_diff, V_INT},
        //{"breakaway_shot_distance",     (void *)&CP_breakaway_shot_distance, V_FLOAT},

        {"formation", (void *)&FP_initial_formation, V_STRING},
        {"formation_when_losing", (void *)&FP_formation_when_losing, V_STRING},
        {"formation_when_losing_lots", (void *)&FP_formation_when_losing_lots, V_STRING},
        {"formation_when_winning", (void *)&FP_formation_when_winning, V_STRING},
        {"formation_when_tied", (void *)&FP_formation_when_tied, V_STRING},

        {"home_change", (void *)&FP_initial_hc_method, V_STRING},
        {"mark_change", (void *)&FP_initial_mc_method, V_STRING},
        {"player_1_pos", (void *)&FP_initial_player_1_pos, V_INT},
        {"player_2_pos", (void *)&FP_initial_player_2_pos, V_INT},
        {"player_3_pos", (void *)&FP_initial_player_3_pos, V_INT},
        {"player_4_pos", (void *)&FP_initial_player_4_pos, V_INT},
        {"player_5_pos", (void *)&FP_initial_player_5_pos, V_INT},
        {"player_6_pos", (void *)&FP_initial_player_6_pos, V_INT},
        {"player_7_pos", (void *)&FP_initial_player_7_pos, V_INT},
        {"player_8_pos", (void *)&FP_initial_player_8_pos, V_INT},
        {"player_9_pos", (void *)&FP_initial_player_9_pos, V_INT},
        {"player_10_pos", (void *)&FP_initial_player_10_pos, V_INT},
        {"player_11_pos", (void *)&FP_initial_player_11_pos, V_INT},
        {"goalie_number", (void *)&FP_goalie_number, V_INT},

        {"version", (void *)&SP_version, V_FLOAT},
        {"team_size", (void *)&SP_team_size, V_INT},
        {"half", (void *)&SP_half, V_INT},
        {"host", (void *)&SP_host, V_STRING},
        {"goal_width", (void *)&SP_goal_width, V_FLOAT},
        {"player_size", (void *)&SP_player_size, V_FLOAT},
        {"player_decay", (void *)&SP_player_decay, V_FLOAT},
        {"player_rand", (void *)&SP_player_rand, V_FLOAT},
        {"player_weight", (void *)&SP_player_weight, V_FLOAT},
        {"player_speed_max", (void *)&SP_player_speed_max, V_FLOAT},
        {"stamina_max", (void *)&SP_stamina_max, V_FLOAT},
        {"stamina_inc_max", (void *)&SP_stamina_inc, V_FLOAT},
        {"recover_dec_thr", (void *)&SP_recover_dec_thr, V_FLOAT},
        {"recover_min", (void *)&SP_recover_min, V_FLOAT},
        {"recover_dec", (void *)&SP_recover_dec, V_FLOAT},
        {"effort_dec_thr", (void *)&SP_effort_dec_thr, V_FLOAT},
        {"effort_min", (void *)&SP_effort_min, V_FLOAT},
        {"effort_dec", (void *)&SP_effort_dec, V_FLOAT},
        {"effort_inc_thr", (void *)&SP_effort_inc_thr, V_FLOAT},
        {"effort_inc", (void *)&SP_effort_inc, V_FLOAT},
        {"ball_size", (void *)&SP_ball_size, V_FLOAT},
        {"ball_decay", (void *)&SP_ball_decay, V_FLOAT},
        {"ball_rand", (void *)&SP_ball_rand, V_FLOAT},
        {"ball_weight", (void *)&SP_ball_weight, V_FLOAT},
        {"ball_speed_max", (void *)&SP_ball_speed_max, V_FLOAT},
        {"dash_power_rate", (void *)&SP_dash_power_rate, V_FLOAT},
        {"kick_power_rate", (void *)&SP_kick_power_rate, V_FLOAT},
        {"kickable_margin", (void *)&SP_kickable_margin, V_FLOAT},
        {"catch_probability", (void *)&SP_catch_prob, V_FLOAT},
        {"catchable_area_l", (void *)&SP_catch_area_l, V_FLOAT},
        {"catchable_area_w", (void *)&SP_catch_area_w, V_FLOAT},
        {"maxpower", (void *)&SP_max_power, V_FLOAT},
        {"minpower", (void *)&SP_min_power, V_FLOAT},
        {"maxmoment", (void *)&SP_max_moment, V_FLOAT},
        {"minmoment", (void *)&SP_min_moment, V_FLOAT},
        {"maxneckang", (void *)&SP_max_neck_angle, V_FLOAT},
        {"minneckang", (void *)&SP_min_neck_angle, V_FLOAT},
        {"maxneckmoment", (void *)&SP_max_neck_moment, V_FLOAT},
        {"minneckmoment", (void *)&SP_min_neck_moment, V_FLOAT},
        {"visible_angle", (void *)&SP_visible_angle, V_FLOAT},
        {"visible_distance", (void *)&SP_visible_dist, V_FLOAT},
        {"audio_cut_dist", (void *)&SP_audio_cut_dist, V_FLOAT},
        {"quantize_step", (void *)&SP_dist_qstep, V_FLOAT},
        {"quantize_step_l", (void *)&SP_land_qstep, V_FLOAT},
        {"ckick_margin", (void *)&SP_ckmargin, V_FLOAT},
        {"wind_dir", (void *)&SP_wind_dir, V_FLOAT},
        {"wind_force", (void *)&SP_wind_force, V_FLOAT},
        {"wind_rand", (void *)&SP_wind_rand, V_FLOAT},
        {"wind_none", (void *)&SP_wind_none, V_ONOFF},
        {"wind_random", (void *)&SP_wind_random, V_ONOFF},
        {"half_time", (void *)&SP_half_time, V_INT},
        {"port", (void *)&SP_port, V_INT},
        {"coach_port", (void *)&SP_coach_port, V_INT},
        {"olcoach_port", (void *)&SP_olcoach_port, V_INT},
        {"simulator_step", (void *)&SP_simulator_step, V_INT},
        {"send_step", (void *)&SP_send_step, V_INT},
        {"recv_step", (void *)&SP_recv_step, V_INT},
        {"say_msg_size", (void *)&SP_say_msg_size, V_INT},
        {"hear_max", (void *)&SP_hear_max, V_INT},
        {"hear_inc", (void *)&SP_hear_inc, V_INT},
        {"hear_decay", (void *)&SP_hear_decay, V_INT},
        {"catch_ban_cycle", (void *)&SP_catch_ban_cycle, V_INT},
        {"coach", (void *)&SP_coach_mode, V_ONOFF},
        {"coach_w_referee", (void *)&SP_coach_w_referee_mode, V_ONOFF},
        {"say_coach_cnt_max", (void *)&SP_say_coach_cnt_max, V_INT},
        {"say_coach_msg_size", (void *)&SP_say_coach_msg_size, V_INT},
        {"send_vi_step", (void *)&SP_send_vi_step, V_INT},

        {"look_step", (void *)&SP_look_step, V_INT},
        {"use_offside", (void *)&SP_use_offside, V_ONOFF},
        {"forbid_kick_off_offside", (void *)&SP_forbid_kickoff_offside, V_ONOFF},
        {"log_file", (void *)&SP_logfile, V_STRING},
        {"record", (void *)&SP_recfile, V_STRING},
        {"record_log", (void *)&SP_rec_log, V_ONOFF},
        {"record_version", (void *)&SP_rec_ver, V_INT},
        {"send_log", (void *)&SP_send_log, V_ONOFF},
        {"replay", (void *)&SP_replay, V_STRING},
        {"verbose", (void *)&SP_verbose, V_ONOFF},
        {"offside_active_area_size", (void *)&SP_offside_area, V_FLOAT},
        {"inertia_moment", (void *)&SP_inertia_moment, V_FLOAT},
        {"sense_body_step", (void *)&SP_sense_body_step, V_INT},
        {"offside_kick_margin", (void *)&SP_offside_kick_margin, V_FLOAT},
        {"record_messages", (void *)&SP_record_messages, V_ONOFF},

        {"\0", NULL, 0}};

    /* skip command name */
    argv++;
    argc--;

    /* first, search option '-file' */

    /* next, analyze command line option */
    int p;

    while (argc)
    {
        if (!strcmp(*argv, "-file"))
        {
            argv += 2;
            argc -= 2;
            continue;
        }

        for (p = 0; opt[p].vptr != NULL; p++)
        {
            if (strcmp(*argv + 1, opt[p].optname))
                continue;

            /* match */
            argv++;
            argc--;

            switch (opt[p].vsize)
            {
            case V_INT:
                *((int *)opt[p].vptr) = atoi(*argv);
                break;

            case V_STRING:
                strcpy((char *)opt[p].vptr, *argv);
                break;

            case V_FLOAT:
                *((float *)opt[p].vptr) = atof(*argv);
                break;

            case V_BOOL:
                *((Bool *)opt[p].vptr) = TRUE;
                argv--;
                argc++;
                break;

            case V_ONOFF:
                if (argc > 0 && (*argv)[0] != '-')
                {
                    *((Bool *)opt[p].vptr) = (!strcmp(*argv, "on")) ? TRUE : FALSE;
                }
                else
                {
                    /* if there's nothing specified, then we set it to true */
                    *((Bool *)opt[p].vptr) = TRUE;
                    argv--;
                    argc++;
                }
                break;
            }

            break;
        }

        if (opt[p].vptr == NULL)
            std::cerr << "Unrecognized Option : " << *argv << std::endl;

        argv++;
        argc--;
    }

    SP_half_time = SP_half_time * 1000 / SP_simulator_step;
    SP_kickable_area = SP_kickable_margin + SP_ball_size + SP_player_size;
}

/* explode the line into argc and argv */
void OptionInfo::GetOptions(char *line)
{
    const int MAXOPT = 100;
    char *argv[MAXOPT];
    int argc = 1; /* executable name */
    char *pc;

    advance_past_space(&line);
    while (*line != 0)
    {
        pc = line;
        get_token(&line);
        argv[argc] = new char[line - pc + 1];
        strncpy(argv[argc], pc, line - pc);
        argv[argc][line - pc] = 0; /* null terminate */
        argc++;
        advance_past_space(&line);
    }

    argv[argc] = NULL;

    GetOptions(argc, argv);

    for (int i = 1; i < argc; i++)
        delete[] argv[i];
}

/* -*- Mode: C++ -*- */

/* MemPlayer.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

PlayerInfo::PlayerInfo()
{
    Initialized = FALSE;
    ServerAlive = FALSE;
    CoachActive = FALSE;
    ViewQuality = VQ_High;
    ViewWidth = VW_Normal;
    NewSight = FALSE;
    NewAction = FALSE;
    FirstActionOpSinceLastSight = FALSE;
    ClockStopped = FALSE;
    StoppedClockMSec = 0;
    LastStartClockTime = Time(-1, 0); /* If a problem, change back to 0,0 */
    SecondLastStartClockTime = LastStartClockTime;
    CurrentTime = 0;
    LastSightTime = LastSoundTime = LastActionOpTime = Time(0, 0);
    PlayModeTime = 0;
    PlayMode = PM_Before_Kick_Off;

    Action = new Command;
    LastAction = new Command;
    RequestResend = FALSE;

    last_dashes = prev_dashes = dashes = 0;
    last_turns = prev_turns = turns = 0;
    last_kicks = prev_kicks = kicks = 0;
    last_says = prev_says = says = 0;
    last_turn_necks = prev_turn_necks = turn_necks = 0;

    TheirTeamName[0] = '\n';

    conf = 0;

    body_ang = 0;
    neck_rel_ang = 0;
}

/*********************************************************************************/

PlayerInfo::~PlayerInfo()
{
    if (CP_save_log)
        fclose(SaveLogFile);

    if (CP_save_sound_log)
        fclose(SaveSoundLogFile);

    if (CP_save_action_log_level > 0)
        fclose(SaveActionLogFile);

    delete Action;
    delete LastAction;
}

/*********************************************************************************/

void PlayerInfo::Initialize()
{
    /* printf("Calling Player Initialize\n"); */

    TheirSide = (MySide == 'l' ? 'r' : 'l');
    MyTeamNameLen = strlen(MyTeamName);

    TestVersion = (VP_test || (MySide == 'l' && VP_test_l) || (MySide == 'r' && VP_test_r)) ? TRUE : FALSE;
    if (TestVersion == TRUE)
        printf("%d : test version\n", MyNumber);
    if (VP_train_DT == TRUE)
        printf("%d : training DT\n", MyNumber);

    MyScore = IP_my_score;
    TheirScore = IP_their_score;

    if (CP_save_log)
    {
        sprintf(SaveLogFileName, "Logfiles/%s%d-%c.log", MyTeamName, (int)MyNumber, MySide);
        SaveLogFile = fopen(SaveLogFileName, "w");
        SaveLogCounter = 0;
    }

    if (CP_save_sound_log)
    {
        sprintf(SaveSoundLogFileName, "Logfiles/%s%d-%c-sounds.log", MyTeamName, (int)MyNumber, MySide);
        SaveSoundLogFile = fopen(SaveSoundLogFileName, "w");
        SaveSoundLogCounter = 0;
    }

    if (CP_save_action_log_level > 0)
    {
        sprintf(SaveActionLogFileName, "Logfiles/%s%d-%c-actions.log", MyTeamName, (int)MyNumber, MySide);
        SaveActionLogFile = fopen(SaveActionLogFileName, "w");
        SaveActionLogCounter = 0;
    }

    TimerInterval = SP_simulator_step / CP_interrupts_per_cycle;

    stamina = SP_stamina_max;
    effort = 1;
    recovery = 1;

    neck_rel_ang = 0;

    RecoveryDecThreshold = SP_stamina_max * SP_recover_dec_thr;
    EffortDecThreshold = SP_stamina_max * SP_effort_dec_thr;
    EffortIncThreshold = SP_stamina_max * SP_effort_inc_thr;

    my_vel_time = my_pos_time = 0;

    LastInterruptTime = 0;
    InterruptsThisCycle = 0;

    Initialized = TRUE;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

#ifndef NO_ACTION_LOG

#define MAX_LOG_LINE 150

/* will be timestamped automatically */
void PlayerInfo::LogAction2(int level, char *str)
{
    if (level <= 0 ||
        level > CP_save_action_log_level)
        return;

    if (!Initialized)
        return; /* the log files hasn't been opened yet! */

    fprintf(Mem->SaveActionLogFile, "%d.%d %s%s\n", CurrentTime.t, CurrentTime.s,
            repeat_char('-', level / 10), str);
    if (Mem->SaveActionLogCounter++ % Mem->CP_save_action_freq == 0)
    {
        fclose(Mem->SaveActionLogFile);
        Mem->SaveActionLogFile = fopen(Mem->SaveActionLogFileName, "a");
    }
}

void PlayerInfo::LogAction3(int level, char *str, char *param)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, param);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction4(int level, char *str, char *param1, char *param2)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, param1, param2);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction3(int level, char *str, char c1)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, c1);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction3(int level, char *str, float f1)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, f1);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction4(int level, char *str, float f1, int d1)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, f1, d1);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction4(int level, char *str, float f1, float f2)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, f1, f2);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction5(int level, char *str, float f1, float f2, float f3)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, f1, f2, f3);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction6(int level, char *str, float f1, float f2, float f3, float f4)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, f1, f2, f3, f4);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction7(int level, char *str, float f1, float f2, float f3, float f4, float f5)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, f1, f2, f3, f4, f5);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction8(int level, char *str, float f1, float f2, float f3, float f4, float f5, float f6)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, f1, f2, f3, f4, f5, f6);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction3(int level, char *str, int d1)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, d1);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction4(int level, char *str, int d1, int d2)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, d1, d2);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction4(int level, char *str, int d1, float f1)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, d1, f1);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction5(int level, char *str, int d1, float f1, float f2)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, d1, f1, f2);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction6(int level, char *str, int d1, float f1, float f2, float f3)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, d1, f1, f2, f3);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction7(int level, char *str, int d1, float f1, float f2, float f3, float f4)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, d1, f1, f2, f3, f4);
    LogAction2(level, outstring);
}

void PlayerInfo::LogAction7(int level, char *str, int d1, int d2, float f1, float f2, float f3)
{
    char outstring[MAX_LOG_LINE];
    sprintf(outstring, str, d1, d2, f1, f2, f3);
    LogAction2(level, outstring);
}

#endif /* ifndef NO_ACTION_LOG */

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

void PlayerInfo::SetPlayMode(Pmode mode)
{
    /* If clock is starting, save the old time */
    if (ClockStopped)
    {
        ClockStopped = FALSE;
        SecondLastStartClockTime = LastStartClockTime;
        LastStartClockTime = LastActionOpTime;
        StoppedClockMSec = 0;

        sanitize_time(CurrentTime);
        sanitize_time(LastSightTime);
        sanitize_time(LastSoundTime);
        sanitize_time(sense_time);
    }

    if ((PlayMode != PM_My_Goalie_Free_Kick && PlayMode != PM_Their_Goalie_Free_Kick) ||
        PlayModeTime != CurrentTime)
    { /* Already set the play mode for this time */
        LastPlayMode = PlayMode;
        PlayMode = mode;
        PlayModeTime = CurrentTime;
    }

    if (mode == PM_Before_Kick_Off || mode == PM_My_Offside_Kick || mode == PM_Their_Offside_Kick)
    {
        if (StoppedClockMSec != 0)
            my_error("StoppedClockMSec should have been reset already");
        ClockStopped = TRUE;
    }

    if (mode == PM_Half_Time || mode == PM_Extended_Time)
        reset_stamina();

    /* in others: */
    //    if ( Mem->QActionTaken )
    //      Mem->CloseRewards();
}

/*********************************************************************************/

void PlayerInfo::sanitize_time(Time &tm)
{
    if (!LastStartClockTime)
        return;

    /* This is to take care of times that were prematurely updated before we knew that
       the clock was about to start again */
    if (tm.t == LastStartClockTime.t && tm.s > LastStartClockTime.s)
    {
        tm = Time(LastStartClockTime.t + (tm.s - LastStartClockTime.s), 0);
    }
}

/*********************************************************************************/

void PlayerInfo::EstimateMyPos()
{
    /* Takes me from previous time to time */
    if (MyConf() && MyVelConf())
        pos += vel;
}

/*********************************************************************************/

void PlayerInfo::EstimateMyVel(Time time)
{
    if (my_vel_time == time && vel_conf == CP_max_conf)
        return;

    /* Takes me from previous time to time */
    if (SensedInfoKnown(time))
    {
        float old_speed = vel.mod();
        if (old_speed)
            vel = vel * GetMySensedSpeed(time) / old_speed; /* don't change the direction */
        else
            vel = Polar2Vector(GetMySensedSpeed(time), MyBodyAng()); /* use my direction */
        vel_conf = CP_max_conf;
    }
    else if (vel_conf < CP_max_conf && SensedInfoKnown(time - 1))
    {
        vel = Polar2Vector(GetMySensedSpeed(time - 1) * SP_player_decay, MyBodyAng());
        vel_conf = CP_conf_decay;
    }
    else if (!MyVelConf())
        return;
    else if (my_vel_time == time - 1)
    {
        vel *= SP_player_decay;
        vel_conf *= CP_conf_decay;
    }
    else if (my_vel_time > time - 10)
    { /* missed up to 10 cycles */
        while (my_vel_time < time && MyVelConf())
        {
            vel *= SP_player_decay;
            vel_conf *= CP_conf_decay;
            ++my_vel_time;
        }
    }
    else
        my_error("Having trouble in vel update -- must have missed at least 10 cycles %.1f %.1f    %f",
                 (float)my_vel_time.t, (float)my_vel_time.s, MyVelConf());

    my_vel_time = time;
}

/*********************************************************************************/

Vector PlayerInfo::NewVelFromDash(Vector old_vel, float dash_power)
{
    float effective_power = MyEffort() * dash_power;
    effective_power *= SP_dash_power_rate;
    Vector new_vel = old_vel + Polar2Vector(effective_power, MyBodyAng());

    if (new_vel.mod() > SP_player_speed_max)
        new_vel *= (SP_player_speed_max / new_vel.mod());

    return new_vel;
}

/*********************************************************************************/

void PlayerInfo::VerifyDash(float *dash_power)
{
    /* Check if recovery going down, or max_speed exceeded */

    float available_power, needed_power = *dash_power;
    if (needed_power < 0)
    {
        needed_power *= -2;
    }
    if (needed_power < 0)
        my_error("power should be positive now");

    float new_stamina = MyStamina() - MyEffort() * needed_power;
    if (new_stamina <= SP_recover_dec_thr * SP_stamina_max && recovery > SP_recover_min)
    {
        /* printf("%d:%d.%d ",MyNumber,CurrentTime.t,CurrentTime.s); */
        /* printf("WARNING: recovery about to go to %.3f\n",recovery - SP_recover_dec); */
        ;
    }
    if (new_stamina <= SP_effort_dec_thr * SP_stamina_max && effort > SP_effort_min)
    {
        /* printf("WARNING: effort about to go to %.2f\n",MyEffort() - SP_effort_dec); */
    }
    if (new_stamina < 0)
    {
        /* printf("%d:%d.%d ",MyNumber,CurrentTime.t,CurrentTime.s); */
        /* printf("WARNING: not enough stamina for dash\n"); */
        available_power = MyStamina() / MyEffort();
        if (*dash_power >= 0)
        {
            *dash_power = available_power;
        }
        else
        {
            *dash_power = -available_power / 2;
        }
    }

    if (NewVelFromDash(MyVel(), *dash_power).mod() > SP_player_speed_max)
    {
        /* printf("%d:%d.%d ",MyNumber,CurrentTime.t,CurrentTime.s); */
        /* printf("WARNING: can't move that fast (assuming vel and dash in same dir)\n"); */
        /* printf("my speed %f   dash_power %f   ",MySpeed(),*dash_power); */
        *dash_power = signf(*dash_power) * (SP_player_speed_max - MySpeed()) / (MyEffort() * SP_dash_power_rate);
        /* printf("new dash_power %f\n",*dash_power); */
    }
}

/*********************************************************************************/

void PlayerInfo::UpdateFromMyAction(Time time)
{
    /* Assume vel and pos are correct for time -- going to time+1 */
    if (!MyConf())
        my_error("Can't update from action if not localized");
    /* But I'm pretty good at estimating... up conf_decay?? */
    if (!NewAction || !(LastActionValid(time)))
        my_error("No action at that time");

    /* AngleDeg delta_ang, expected_delta; */
    switch (LastActionType())
    {
    case CMD_turn:
        if (my_pos_time > time)
            break;
        /* be careful not to estimate in a turn that's already been seen --
           server updates turns instantaneously */
        /* THIS SHOULDN'T HAPPEN ANYMORE */
        /*     delta_ang = GetNormalizeAngleDeg(ang - my_last_ang); */
        /*     expected_delta = LastActionAngle()/(1.0 + SP_inertia_moment * MySpeed()); */

        /* only if the change is closer to 0 than to the expected change */
        /*     if ( fabs(delta_ang) < fabs(delta_ang-expected_delta) ){ */
        /*        body_ang += expected_delta; */
        body_ang += LastActionAngle() / (1.0 + SP_inertia_moment * MySpeed());
        NormalizeAngleDeg(&body_ang);
        /*     } */
        /*     else */
        /*       my_error("Turns should NOT happen instantaneously anymore"); */
        break;
    case CMD_dash:
        if (my_vel_time > time)
            break;
        vel = NewVelFromDash(vel, LastActionPower());
        break;
    default:;
    }
}

/*********************************************************************************/

void PlayerInfo::update_self_estimate(Time time)
{
    update_self_neck_rel_ang(time);

    if (!MyConf())
    {
        vel_conf = 0; /* If don't know my position, can't know my velocity */
        return;
    }

    if (CP_use_new_position_based_vel)
    {
        if (my_pos_time == time)
        { /* just vel */
            if (my_vel_time == time)
                return;
            if (NewAction && LastActionValid(my_vel_time))
                UpdateFromMyAction(my_vel_time);

            EstimateMyVel(time);
        }
    }
    else
    {
        if (my_pos_time == time)
        { /* just vel */
            if (my_vel_time == time)
                my_error("my pos and vel already updated\n");
            if (NewAction && LastActionValid(my_vel_time))
                UpdateFromMyAction(my_vel_time);

            EstimateMyVel(time);
        }
    }

    while (my_pos_time < time)
    {
        if (NewAction && LastActionValid(my_pos_time))
            UpdateFromMyAction(my_pos_time);

        ++my_pos_time;

        EstimateMyPos();
        EstimateMyVel(time);

        conf *= CP_conf_decay;
    }
}

/*********************************************************************************/

void PlayerInfo::update_self_neck_rel_ang(Time time)
{

    if (SensedInfoKnown(time))
        SetMyNeckRelAng(GetMySensedNeckAngle(time));
    else if (SensedInfoKnown(time - 1))
    {
        /* Bring it up to date from the action */
        AngleDeg neck_ang = GetMySensedNeckAngle(time - 1);
        if (TurnNeck.valid(time - 1))
        {
            neck_ang += TurnNeck.angle;
            if (neck_ang < SP_min_neck_angle)
                neck_ang = SP_min_neck_angle;
            if (neck_ang > SP_max_neck_angle)
                neck_ang = SP_max_neck_angle;
        }
        SetMyNeckRelAng(neck_ang);
    }
    else
        /* could write an "estimate_neck_angle" that updates from the last known time */
        /* could also assume neck unchanged */
        ; /*my_error("Don't know neck angle at time %d.%d or %d.%d",
              time.t,time.s,(time-1).t,(time-1).s);*/
}

/*********************************************************************************/

void PlayerInfo::update_stamina(Time time)
{
    if (NewAction && LastActionType() == CMD_dash)
        stamina -= (LastActionPower() > 0) ? LastActionPower() : (-2.0 * LastActionPower());

    if (stamina < 0)
        stamina = 0;

    if (stamina <= SP_recover_dec_thr * SP_stamina_max && recovery > SP_recover_min)
    {
        recovery -= SP_recover_dec;
    }

    if (SensedInfoKnown(time))
    {
        stamina = GetMySensedStamina(time);
        effort = GetMySensedEffort(time);
    }
    else
    {
        if (stamina <= SP_effort_dec_thr * SP_stamina_max && effort > SP_effort_min)
            effort -= SP_effort_dec;
        if (stamina >= SP_effort_inc_thr * SP_stamina_max && effort < 1.0)
        {
            effort += SP_effort_inc;
            if (effort > 1.0)
                effort = 1.0;
        }
        stamina += recovery * SP_stamina_inc;
        if (stamina > SP_stamina_max)
            stamina = SP_stamina_max;
    }
}

/*********************************************************************************/

void PlayerInfo::reset_stamina()
{
    stamina = SP_stamina_max;
    effort = recovery = 1.0;
}

/*********************************************************************************/

Time PlayerInfo::update_time(int time)
{
    LastTime = CurrentTime;

    if (ClockStopped)
    {
        if (CurrentTime.t != time)
        {
            if (CurrentTime.t == time - 1) /* Sometimes happens in offsides mode */
                CurrentTime = Time(time, 0);
            else
                my_error("server time should be the same %d %d %d", CurrentTime.t, CurrentTime.s, time);
        }
        else
            CurrentTime.s = StoppedClockMSec / SP_simulator_step;
    }
    else if (LastStartClockTime.t == time)
        CurrentTime = LastStartClockTime;
    else
        CurrentTime = Time(time, 0);

    return CurrentTime;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

Bool PlayerInfo::SightPredictedEarlyThisCycle()
{
    if (InterruptsThisCycle > CP_interrupts_per_cycle / 2)
        /* already past the beginning of the cycle */
        return FALSE;

    /* Number of full cycles since last sight * simulator_step */
    if (MySightInterval() - ((CurrentTime - LastSightTime) - 1) * SP_simulator_step <= SP_simulator_step / 2)
        return TRUE;

    return FALSE;
}

/*********************************************************************************/

Bool PlayerInfo::GotSightFromCurrentPosition()
{
    if (FirstActionOpSinceLastSight &&
        /* sight from this time or didn't change view angle last time step */
        /* could replace valids with MyNeckGlobalAng() == my_last_neck_global_ang)) */
        /* but when the angle's estimated, it might be off by up to 10 degrees */
        (LastSightTime == CurrentTime ||
         (!LastAction->valid(CurrentTime - 1) && !TurnNeck.valid(CurrentTime - 1))))
        return TRUE;

    return FALSE;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

AngleDeg PlayerInfo::MyViewAngle(Time time)
{
    AngleDeg view_angle = SP_visible_angle;
    Vwidth width;

    if (time < ViewWidthTime)
        width = LastViewWidth;
    else
        width = ViewWidth;

    if (width == VW_Narrow)
        view_angle /= 2;
    else if (width == VW_Wide)
        view_angle *= 2;

    return view_angle / 2;
}

/*********************************************************************************/

Bool PlayerInfo::InViewAngle(Time time, AngleDeg ang, float buffer)
{
    if (fabs(ang) < MyViewAngle(time) - buffer)
        return TRUE;
    return FALSE;
}

/*********************************************************************************/

int PlayerInfo::MySightInterval()
{
    int interval = SP_send_step;

    if (ViewWidth == VW_Narrow)
        interval /= 2;
    else if (ViewWidth == VW_Wide)
        interval *= 2;

    if (ViewQuality == VQ_Low)
        interval /= 2;

    return interval;
}

/*********************************************************************************/

int PlayerInfo::PredictedNextSightInterval()
{
    int interval = MySightInterval();
    if (interval < SP_simulator_step) /* 37 or 75 */
        return 1;
    if (interval == 3 * SP_simulator_step) /* 300 */
        return 3;
    if (interval == 1.5 * SP_simulator_step) /* 150 */
        return (LastSightInterval <= 1 ? 2 : 1);

    my_error("Sight interval should be 37, 75, 150, or 300: %d", MySightInterval());
    return 0;
}

/*********************************************************************************/

void PlayerInfo::SetMySensedInfo(float st, float e, float sp, float ha, int k, int d, int tu, int sa, int tn, Time ti)
{
    if (sense_time == ti)
        return;

    prev_sense_time = sense_time;
    sense_time = ti;

    prev_stamina = last_stamina;
    last_stamina = st;
    prev_effort = last_effort;
    last_effort = e;
    prev_speed = last_speed;
    last_speed = sp;
    prev_neck_rel_ang = last_neck_rel_ang;
    last_neck_rel_ang = ha;

    //  neck_rel_ang    = ha;  /** Want to do this here??? No! **/

    prev_kicks = last_kicks;
    last_kicks = k;
    if (last_kicks != kicks)
    {
        if (!ClockStopped)
            my_error("Server missed a kick at time %d (%d %d)", prev_sense_time.t, last_kicks, kicks);
        LastAction->type = CMD_none;
        kicks = last_kicks;
        Mem->GetBall()->forget_past_kick(LastAction->time);
        /* RequestResend = TRUE;
           ResendType    = CMD_kick;
           ResendTime    = LastActionTime(); */
    }

    prev_dashes = last_dashes;
    last_dashes = d;
    if (last_dashes != dashes)
    {
        if (!ClockStopped)
            my_error("Server missed a dash at time %d (%d %d)", prev_sense_time.t, last_dashes, dashes);
        LastAction->type = CMD_none;
        dashes = last_dashes;
        /* RequestResend = TRUE;
           ResendType   = CMD_dash;
           ResendTime    = LastActionTime(); */
    }

    prev_turns = last_turns;
    last_turns = tu;
    if (last_turns != turns)
    {
        if (!ClockStopped)
            my_error("Server missed a turn at time %d (%d %d)", prev_sense_time.t, last_turns, turns);
        LastAction->type = CMD_none;
        turns = last_turns;
        /* RequestResend = TRUE;
           ResendType   = CMD_turn;
           ResendTime    = LastActionTime(); */
    }

    prev_turn_necks = last_turn_necks;
    last_turn_necks = tn;
    if (last_turn_necks != turn_necks)
    {
        if (!ClockStopped)
            my_error("Server missed a turn_neck at time %d (%d %d)", prev_sense_time.t, last_turn_necks, turn_necks);
        TurnNeck.type = CMD_none;
        turn_necks = last_turn_necks;
        /* RequestResend = TRUE;
           ResendType   = CMD_turn;
           ResendTime    = LastActionTime(); */
    }

    prev_says = last_says;
    last_says = sa;
    if (last_says != says)
    {
        says = last_says;
    }
}

/*********************************************************************************/

float PlayerInfo::GetMySensedSpeed(Time time)
{

    if (time == sense_time)
        return last_speed;
    if (time == prev_sense_time)
        return prev_speed;

    my_error("Don't know my speed at time %d", time.t);
    return 0;
}

/*********************************************************************************/

float PlayerInfo::GetMySensedStamina(Time time)
{

    if (time == sense_time)
        return last_stamina;
    if (time == prev_sense_time)
        return prev_stamina;

    my_error("Don't know my stamina at time %d", time.t);
    return 0;
}

/*********************************************************************************/

float PlayerInfo::GetMySensedEffort(Time time)
{

    if (time == sense_time)
        return last_effort;
    if (time == prev_sense_time)
        return prev_effort;

    my_error("Don't know my effort at time %d", time.t);
    return 0;
}

/*********************************************************************************/

float PlayerInfo::GetMySensedNeckAngle(Time time)
{

    if (time == sense_time)
        return last_neck_rel_ang;
    if (time == prev_sense_time)
        return prev_neck_rel_ang;

    my_error("Don't know my neck angle at time %d", time.t);
    return 0;
}

/*********************************************************************************/

int PlayerInfo::GetMySensedKicks(Time time)
{

    if (time == sense_time)
        return last_kicks;
    if (time == prev_sense_time)
        return prev_kicks;

    my_error("Don't know my kicks at time %d", time.t);
    return 0;
}

/*********************************************************************************/

int PlayerInfo::GetMySensedDashes(Time time)
{

    if (time == sense_time)
        return last_dashes;
    if (time == prev_sense_time)
        return prev_dashes;

    my_error("Don't know my dashes at time %d", time.t);
    return 0;
}

/*********************************************************************************/

int PlayerInfo::GetMySensedTurns(Time time)
{

    if (time == sense_time)
        return last_turns;
    if (time == prev_sense_time)
        return prev_turns;

    my_error("Don't know my turns at time %d", time.t);
    return 0;
}

/*********************************************************************************/

int PlayerInfo::GetMySensedSays(Time time)
{

    if (time == sense_time)
        return last_says;
    if (time == prev_sense_time)
        return prev_says;

    my_error("Don't know my says at time %d", time.t);
    return 0;
}

/*********************************************************************************/

int PlayerInfo::GetMySensedTurnNecks(Time time)
{

    if (time == sense_time)
        return last_turn_necks;
    if (time == prev_sense_time)
        return prev_turn_necks;

    my_error("Don't know my turn_necks at time %d", time.t);
    return 0;
}

/*********************************************************************************/

float PlayerInfo::CorrectDashPowerForStamina(float dash_power, float stamina, float, float)
{
    float new_power;
    if (dash_power >= 0)
    {
        new_power = Min(dash_power, stamina - (EffortDecThreshold + CP_tired_buffer));
        if (new_power < 0)
            new_power = 0;
    }
    else
    {
        new_power = Min(-dash_power, stamina - (EffortDecThreshold + CP_tired_buffer) / 2.0);
        if (new_power < 0)
            new_power = 0;

        new_power = -new_power;
    }

    return new_power;
}

/*********************************************************************************/

Bool PlayerInfo::CanFaceAngleFromNeckWithNeck(AngleDeg ang)
{
    AngleDeg total_ang = MyNeckRelAng() + ang;
    NormalizeAngleDeg(&total_ang);
    if (total_ang > SP_min_neck_angle && total_ang < SP_max_neck_angle)
        return TRUE;
    return FALSE;
}

/*********************************************************************************/

Bool PlayerInfo::CanFaceAngleFromBodyWithNeck(AngleDeg ang)
{
    NormalizeAngleDeg(&ang);
    if (ang > SP_min_neck_angle && ang < SP_max_neck_angle)
        return TRUE;
    return FALSE;
}

/*********************************************************************************/

Bool PlayerInfo::CanSeeAngleFromNeckWithNeck(AngleDeg ang)
{
    AngleDeg total_ang = MyNeckRelAng() + ang;
    NormalizeAngleDeg(&total_ang);
    if (total_ang > SP_min_neck_angle - MyViewAngle() &&
        total_ang < SP_max_neck_angle + MyViewAngle())
        return TRUE;
    return FALSE;
}

/*********************************************************************************/

Bool PlayerInfo::CanSeeAngleFromBodyWithNeck(AngleDeg ang)
{
    if (ang > 180 || ang < -180)
    {
        my_error("Passing unnormalized angle to CanSeeAngleFromBodyWithNeck: %.1f", ang);
        NormalizeAngleDeg(&ang);
    }
    if (ang > SP_min_neck_angle - MyViewAngle() &&
        ang < SP_max_neck_angle + MyViewAngle())
        return TRUE;
    return FALSE;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
void PlayerInfo::UpdatePredictedStaminaWithDash(float *pStamina, float *pEffort,
                                                float *pRecovery, float dash_power)
{
    if (dash_power > 0)
        *pStamina -= dash_power;
    else
        *pStamina -= 2 * dash_power;
    if (*pStamina < 0)
        *pStamina = 0;

    if (*pStamina <= SP_recover_dec_thr * SP_stamina_max && *pRecovery > SP_recover_min)
    {
        *pRecovery -= SP_recover_dec;
    }

    if (*pStamina <= SP_effort_dec_thr * SP_stamina_max && *pEffort > SP_effort_min)
        *pEffort -= SP_effort_dec;
    if (*pStamina >= SP_effort_inc_thr * SP_stamina_max && *pEffort < 1.0)
    {
        *pEffort += SP_effort_inc;
        if (*pEffort > 1.0)
            *pEffort = 1.0;
    }
    *pStamina += *pRecovery * SP_stamina_inc;
    if (*pStamina > SP_stamina_max)
        *pStamina = SP_stamina_max;
}

/*********************************************************************************/

Vector PlayerInfo::MyPredictedPositionAtMaxSpeed(int steps)
{
    if (!MyConf())
        my_error("Can't estimate future if don't know present (max speed)");

    Vector new_position = MyPos();
    Vector max_velocity = Polar2Vector(SP_player_speed_max, MyBodyAng());
    for (int i = 0; i < steps; i++)
    {
        new_position += max_velocity;
    }
    return new_position;
}

/*********************************************************************************/

Vector PlayerInfo::MyPredictedPositionWithTurn(float turn_ang,
                                               int steps, float dash_power,
                                               bool with_turn,
                                               int idle_cycles)
{
    if (!MyConf())
        my_error("Can't estimate future if don't know present");

    float curr_turn_ang = GetNormalizeAngleDeg(turn_ang);
    float corrected_dash_power = dash_power;
    float effective_power;
    float predicted_stamina = MyStamina();
    float predicted_effort = MyEffort();
    float predicted_recovery = MyRecovery();
    float myang = MyBodyAng();
    Vector position = MyPos();
    Vector velocity;
    if (!MyVelConf())
        velocity = 0;
    else
        velocity = MyVel();
    /* debug code
    cout << "steps: " << steps << "\tpow: " << dash_power << "\tmyang: " << myang
         << "\tposition: " << position << "\tvel: " << velocity
         << "\tturn?: " << turn_first << "\tturn_ang: " << turn_angle
         << "\tstam: " << predicted_stamina << "\teff: " << predicted_effort
         << "\trec: " << predicted_recovery << endl; */

    for (int i = 0; i < steps; i++)
    {
        corrected_dash_power = CorrectDashPowerForStamina(dash_power, predicted_stamina);
        /* cout << " in func: i=" << i << "\tpos" << position << endl; */
        if (i < idle_cycles)
        {
            /* do nothing, we're idling! */
            effective_power = 0;
        }
        else if (with_turn &&
                 (i == 0 || curr_turn_ang != 0.0))
        {
            float this_turn = MinMax(-EffectiveTurn(SP_max_moment, velocity.mod()),
                                     curr_turn_ang,
                                     EffectiveTurn(SP_max_moment, velocity.mod()));
            myang += this_turn;
            curr_turn_ang -= this_turn;
            effective_power = 0;
        }
        else if (fabs(corrected_dash_power) > predicted_stamina)
            effective_power = Sign(corrected_dash_power) * predicted_stamina;
        else
            effective_power = corrected_dash_power;

        effective_power *= predicted_effort;
        effective_power *= SP_dash_power_rate;
        velocity += Polar2Vector(effective_power, myang);
        /* cout << " in func: i=" << i << "\tvel" << velocity << endl; */

        if (velocity.mod() > SP_player_speed_max)
            velocity *= (SP_player_speed_max / velocity.mod());

        position += velocity;
        velocity *= SP_player_decay;

        UpdatePredictedStaminaWithDash(&predicted_stamina, &predicted_effort,
                                       &predicted_recovery, corrected_dash_power);

        /*
        predicted_stamina -= predicted_effort * fabs(corrected_dash_power);
        if (predicted_stamina < 0) predicted_stamina = 0;

        if ( predicted_stamina <= SP_recover_dec_thr * SP_stamina_max && predicted_recovery > SP_recover_min ) {
          predicted_recovery -= SP_recover_dec;
        }

        if ( predicted_stamina <= SP_effort_dec_thr * SP_stamina_max && predicted_effort > SP_effort_min )
          predicted_effort -= SP_effort_dec;
        if (predicted_stamina >= SP_effort_inc_thr * SP_stamina_max && predicted_effort < 1.0){
          predicted_effort += SP_effort_inc;
          if ( predicted_effort > 1.0 )
        predicted_effort = 1.0;
        }
        predicted_stamina += predicted_recovery * SP_stamina_inc;
        if ( predicted_stamina > SP_stamina_max )
          predicted_stamina = SP_stamina_max;
      */
    }
    /* cout << "returning " << position << endl; */
    return position;
}

/*********************************************************************************/

Vector PlayerInfo::MyPredictedPositionWithQueuedActions()
{
    /* Only goes one step in the future so far (other function assumes repeated dashes) */
    if (Action->valid() && Action->type == CMD_dash)
        return MyPredictedPosition(1, Action->power);
    else
        return MyPredictedPosition();
}

/*********************************************************************************/

AngleDeg PlayerInfo::MyPredictedBodyAngleWithQueuedActions()
{
    /* Only goes one step in the future so far (other function assumes repeated dashes) */
    if (Action->valid() && Action->type == CMD_turn)
        return GetNormalizeAngleDeg(MyBodyAng() + EffectiveTurn(Action->angle));
    else
        return MyBodyAng();
}

/********************************************************************************/

AngleDeg PlayerInfo::PredictedPointRelAngFromBodyWithQueuedActions(Vector point)
{
    Vector pred_my_pos = MyPredictedPositionWithQueuedActions();
    AngleDeg pred_my_body_ang = MyPredictedBodyAngleWithQueuedActions();
    Vector pred_relative_point_pos = point - pred_my_pos;
    AngleDeg target_abs_ang = pred_relative_point_pos.dir();
    AngleDeg target_rel_ang = target_abs_ang - pred_my_body_ang;
    NormalizeAngleDeg(&target_rel_ang);

    return target_rel_ang;
}

/*********************************************************************************/

int PlayerInfo::PredictedCyclesToPoint(Vector pt, float dash_power)
{
    float corrected_dash_power = dash_power;
    float effective_power;
    float predicted_stamina = MyStamina();
    float predicted_effort = MyEffort();
    float predicted_recovery = MyRecovery();
    float myang = MyBodyAng();
    Vector position = MyPos();
    Vector velocity;
    if (!MyVelConf())
        velocity = 0;
    else
        velocity = MyVel();

    for (int i = 0; TRUE; i++)
    {
        if (position.dist(pt) <= CP_at_point_buffer)
            return i;

        /* decide if we should turn */
        float targ_ang = (pt - position).dir() - myang;
        if (fabs(GetNormalizeAngleDeg(targ_ang)) > CP_max_go_to_point_angle_err)
        {
            /* turning */
            float this_turn = MinMax(-EffectiveTurn(SP_max_moment, velocity.mod()),
                                     targ_ang,
                                     EffectiveTurn(SP_max_moment, velocity.mod()));
            myang += this_turn;
            corrected_dash_power = 0; // so that stamina is updated correctly
        }
        else
        {
            /* dashing */
            corrected_dash_power = CorrectDashPowerForStamina(dash_power, predicted_stamina);
            if (fabs(corrected_dash_power) > predicted_stamina)
                effective_power = Sign(corrected_dash_power) * predicted_stamina;
            else
                effective_power = corrected_dash_power;

            effective_power *= predicted_effort;
            effective_power *= SP_dash_power_rate;
            velocity += Polar2Vector(effective_power, myang);
        }

        if (velocity.mod() > SP_player_speed_max)
            velocity *= (SP_player_speed_max / velocity.mod());

        position += velocity;
        velocity *= SP_player_decay;

        UpdatePredictedStaminaWithDash(&predicted_stamina, &predicted_effort,
                                       &predicted_recovery, corrected_dash_power);
    }
}

/*********************************************************************************/
int PlayerInfo::NumTurnsToAngle(float targ_body_ang, float curr_body_ang, float curr_speed)
{
    int steps;

    NormalizeAngleDeg(&targ_body_ang);
    NormalizeAngleDeg(&curr_body_ang);

    for (steps = 0;
         fabs(targ_body_ang - curr_body_ang) > CP_max_go_to_point_angle_err;
         steps++)
    {
        AngleDeg this_turn = targ_body_ang - curr_body_ang;
        NormalizeAngleDeg(&this_turn);
        this_turn = signf(this_turn) * Min(fabs(this_turn), MaxEffectiveTurn(curr_speed));
        Mem->LogAction5(210, "NumTurnsToAngle: curr: %.1f  targ: %.1f  turn: %.1f",
                        curr_body_ang, targ_body_ang, this_turn);
        curr_body_ang += this_turn;
        NormalizeAngleDeg(&curr_body_ang);
        curr_speed *= SP_player_decay;
    }

    return steps;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/
/* returns whether we're trying to send a command too soon after the previous command.
   Helps to keep us from missing commands */
/* Not needed after server 5.23
Bool PlayerInfo::TooSoonForAnotherSend()
{
  struct timeval tv_new;

  gettimeofday(&tv_new, NULL); // no time zone info;
  int usec_diff = (tv_new.tv_sec - real_time_of_last_send.tv_sec) * 1000000 +
    ((signed)tv_new.tv_usec - (signed)real_time_of_last_send.tv_usec);
  return (usec_diff < Mem->SP_recv_step*1000 / Mem->CP_send_ban_recv_step_factor)
    ? TRUE : FALSE;
}
*/

/* -*- Mode: C++ -*- */

/* MemPosition.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        Object Class                                          */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

float Object::get_dist()
{
    if (pos_valid() && Mem->MyConf())
        return Mem->DistanceTo(gpos); /* From global   */

    /* dump_core("dump"); */
    my_error("Can't get_dist %d", type);
    return 0;

    if (pdtime < Mem->CurrentTime - 1) /* In sight      */
        my_error("Can't get_dist %d", type);

    printf("Using in-sight estimate of dist\n");
    return dist; /* From relative */
}

/********************************************************************************/

AngleDeg Object::get_ang_from_body()
{
    if (pos_valid() || Mem->MyConf())
        return Mem->AngleToFromBody(gpos); /* From global   */

    my_error("Can't get_ang %d", type);
    return 0;

    if (patime < Mem->CurrentTime - 1)
        my_error("Can't get_ang %d", type);

    /* In sight */
    printf("Using in-sight estimate of ang\n");
    return ang_from_neck + Mem->MyNeckRelAng(); /* From relative */
}

/********************************************************************************/

AngleDeg Object::get_ang_from_neck()
{
    if (pos_valid() || Mem->MyConf())
        return Mem->AngleToFromNeck(gpos); /* From global   */

    my_error("Can't get_ang %d", type);
    return 0;

    if (patime < Mem->CurrentTime - 1)
        my_error("Can't get_ang %d", type);

    /* In sight */
    printf("Using in-sight estimate of ang\n");
    return ang_from_neck; /* From relative */
}

/********************************************************************************/
/* we have some problems keeping these consitent I think,
   (especially during the update process, which has the positions in a state)
   that is nto equal to the current time
   therefore, we'll just recompute them all the time */

Vector Object::get_rel_to_body_pos() /* relative */
{
    return (gpos - Mem->MyPos()).rotate(-Mem->MyBodyAng());
    /* see the somment above
    if ( rbtime != Mem->CurrentTime ){
      rbpos = (gpos - Mem->MyPos()).rotate(-Mem->MyBodyAng());
      rbtime = Mem->CurrentTime;
    }
    return rbpos;
    */
}

/********************************************************************************/

Vector Object::get_rel_to_neck_pos() /* relative */
{
    return (gpos - Mem->MyPos()).rotate(-Mem->MyNeckGlobalAng());
    /* see the somment above
    if ( rntime != Mem->CurrentTime ){
      rnpos = (gpos - Mem->MyPos()).rotate(-Mem->MyNeckGlobalAng());
      rntime = Mem->CurrentTime;
    }
    return rnpos;
    */
}

/********************************************************************************/

void Object::set_polar_from_neck(float d, float a, Time time)
{
    dist = d;
    ang_from_neck = a;
    pdtime = patime = time;
    seen = TRUE;
    seen_time = time;
}

/********************************************************************************/

void Object::set_angle_from_neck(AngleDeg a, Time time)
{
    ang_from_neck = a;
    patime = time;
    seen = TRUE;
    seen_time = time;
}

/********************************************************************************/

void Object::set_chinfo(float dist, float dir, Time time)
{
    distch = dist;
    dirch = dir;
    chtime = time;
    seen_moving = TRUE;
}

/********************************************************************************/

void Object::update()
{
    seen = seen_moving = FALSE;
}

/********************************************************************************/

void Object::reset()
{
    chtime = pdtime = patime = 0;
    seen = FALSE;
    seen_time = -3;
}

/********************************************************************************/

void Object::clear_seen()
{
    seen = FALSE;
}

/********************************************************************************/

void Object::sanitize_times()
{
    Mem->sanitize_time(seen_time);
    Mem->sanitize_time(pdtime);
    Mem->sanitize_time(patime);
    Mem->sanitize_time(chtime);
}

/********************************************************************************/

Bool Object::in_view_range(AngleDeg view_ang, float angle_buffer, float distance_buffer)
{
    if (!pos_valid() || !Mem->MyConf())
        return FALSE;
    if (distance_buffer != 0)
        my_error("Shouldn't have valid distance_buffer at this point");

    Vector tmp_pos = get_rel_to_neck_pos();
    /* angle_buffer is the maximum distance that an object should be off to be forgotten
       We slide the object back before taking the angle in order to accomplish this */
    float slide_dist = angle_buffer / Sin(view_ang);
    // tmp_pos -= Vector(angle_buffer,0);
    tmp_pos -= Vector(slide_dist, 0);

    if (fabs(tmp_pos.dir()) <= view_ang)
    {
        Mem->LogAction7(200, "Object::in_view_range, tmp_pos: (%.2f, %.2f) %.2f %.2f %.2f",
                        tmp_pos.x, tmp_pos.y, tmp_pos.dir(), view_ang, slide_dist);
        return TRUE;
    }

    /* if ( fabs(get_ang()) <= Mem->MyViewAngle() - angle_buffer ) return TRUE; */

    return FALSE;
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        StationaryObject Class                                */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void StationaryObject::Initialize(MarkerType m, Vector pos, float max, float min_valid, Bool rotate)
{
    gpos = rotate ? -pos : pos; /* If on the right side of the field, flip all coords */

    max_conf = max;
    min_valid_conf = min_valid;

    gconf = max_conf;

    type = OBJ_Marker;
    object_id = m;
}

/********************************************************************************/

void StationaryObject::Initialize(SideLine l, Vector pos, float max, float min_valid, Bool rotate)
{
    gpos = rotate ? -pos : pos; /* If on the right side of the field, flip all coords */

    max_conf = max;
    min_valid_conf = min_valid;

    gconf = max_conf;

    type = OBJ_Line;
    object_id = l;
}

/********************************************************************************/

Vector StationaryObject::get_my_pos(AngleDeg my_neck_global_ang)
{
    if (type != OBJ_Marker)
        my_error("Need to get pos with a marker\n");

    AngleDeg abs_ang = ang_from_neck + my_neck_global_ang;
    NormalizeAngleDeg(&abs_ang);
    Vector rpos = Polar2Vector(dist, abs_ang);

    return gpos - rpos;
}

/********************************************************************************/

AngleDeg StationaryObject::get_my_neck_global_ang()
{

    if (type != OBJ_Line)
    {
        my_error("Need to get angle with a line\n");
    }

    AngleDeg line_ang = gpos.dir();
    AngleDeg my_neck_ang;

    my_neck_ang = (ang_from_neck < 0) ? line_ang - (90 + ang_from_neck) : line_ang + (90 - ang_from_neck);

    NormalizeAngleDeg(&my_neck_ang);

    return my_neck_ang;
}

/********************************************************************************/

Vector StationaryObject::get_my_vel(AngleDeg my_neck_global_ang)
{
    my_error("Shouldn't be estimating velocity from markers -- use sense_body");
    if (type != OBJ_Marker)
        my_error("Need to get vel with a marker\n");

    Vector rvel = Vector(distch, dist * Tan(dirch)).rotate(ang_from_neck + my_neck_global_ang);
    return -rvel; /* Assume the object's not moving */
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        MobileObject Class                                    */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void MobileObject::Initialize(ObjType t, float max, float min_valid, float decay, float motion, float speed)
{
    max_conf = max;
    min_valid_conf = min_valid;
    conf_decay = decay;
    max_speed = speed;
    motion_decay = motion;
    type = t;
    reset();
}

/********************************************************************************/

AngleDeg MobileObject::get_rel_to_body_heading()
{
    float h = get_abs_heading() - Mem->MyBodyAng();
    NormalizeAngleDeg(&h);
    return h;
}

/********************************************************************************/

AngleDeg MobileObject::get_rel_to_neck_heading()
{
    float h = get_abs_heading() - Mem->MyNeckGlobalAng();
    NormalizeAngleDeg(&h);
    return h;
}

/********************************************************************************/

Vector MobileObject::get_rel_to_body_vel() /* relative */
{
    if (rbvtime != Mem->CurrentTime)
    {
        /* rvel = (gvel - Mem->MyVel()).rotate(-Mem->MyAng()); */ /* TRUE RELATIVE VEL            */
        rbvel = gvel.rotate(-Mem->MyBodyAng());                   /* ABSOLUTE SPEED, RELATIVE DIR */
        rbvtime = Mem->CurrentTime;
    }
    return rbvel;
}

/********************************************************************************/

Vector MobileObject::get_rel_to_neck_vel() /* relative */
{
    if (rnvtime != Mem->CurrentTime)
    {
        /* rvel = (gvel - Mem->MyVel()).rotate(-Mem->MyAng()); */ /* TRUE RELATIVE VEL            */
        rnvel = gvel.rotate(-Mem->MyNeckGlobalAng());             /* ABSOLUTE SPEED, RELATIVE DIR */
        rnvtime = Mem->CurrentTime;
    }
    return rnvel;
}

/********************************************************************************/

void MobileObject::set_heard_info(float x, float y, float conf, float dist, Time time)
{
    /* Don't overwrite better infor from the same cycle */
    if (heard &&
        (conf < hpconf || (conf == hpconf && hpdist > dist)))
        return;

    hpos = Vector(x, y);
    hpconf = conf;
    hpdist = dist;
    hptime = time;
    heard = TRUE;
}

/********************************************************************************/

void MobileObject::set_heard_info(float x, float y, float pconf, float dx, float dy, float vconf, float dist, Time time)
{
    /* Don't overwrite better infor from the same cycle      */
    /* If I have better motion info, I have better pos. info */
    if (heard_moving &&
        (vconf < hvconf || (vconf == hvconf && hvdist > dist)))
        return;

    set_heard_info(x, y, pconf, dist, time);
    hvel = Vector(dx, dy);
    hvconf = vconf;
    hvdist = dist;
    hvtime = time;
    heard_moving = TRUE;

    float speed = hvel.mod();
    if (speed > max_speed)
    {
        if (speed > max_speed + .1)
            my_error("object CAN'T be moving that fast %.0f %.2f", (float)type, speed);
        hvel *= max_speed / speed;
    }
}

/********************************************************************************/

void MobileObject::estimate_pos(Time time)
{
    if (gtime >= time)
        my_error("pos already updated");

    if (pos_valid() && vel_valid())
        gpos += gvel;
}

/********************************************************************************/

void MobileObject::estimate_vel(Time time)
{
    if (gvtime >= time)
    {
        my_error("vel already updated");
    }

    if (vel_valid())
        gvel *= motion_decay;
    else
        my_error("Shouldn't be updating invalid vel");
}

/********************************************************************************/

void MobileObject::update(Time time)
{
    if (Mem->NewSight && (seen || seen_moving))
        update_seen(Mem->LastSightTime);

    /* change_view happens instantaneously, so last sight could've been from either view angle */
    AngleDeg forget_view_ang = Min(Mem->MyViewAngle(time), Mem->MyViewAngle(time - 1));
    float ang_buff = (type == OBJ_Ball ? Mem->CP_ball_forget_angle_buf : Mem->CP_player_forget_angle_buf);
    float dist_buff = (type == OBJ_Ball ? Mem->CP_ball_forget_dist_buf : Mem->CP_player_forget_dist_buf);

    /* If sight was from time-1, and should be in view at this point but isn't, reset */
    if (Mem->NewSight && Mem->LastSightTime == time - 1 &&
        !seen && in_view_range(forget_view_ang, ang_buff, dist_buff))
        forget();

    if (gvtime < time || gtime < time)
        update_estimate(time);

    /* If sight was from time, and should be in view at this point but isn't, reset */
    if (Mem->NewSight && Mem->LastSightTime == time &&
        !seen && in_view_range(forget_view_ang, ang_buff, dist_buff))
        forget();

    if (heard || heard_moving)
        update_heard(time);

    Object::update();
    heard = heard_moving = FALSE;
}

/********************************************************************************/

void MobileObject::update_seen(Time time)
{
    if (!Mem->MyConf())
        return;

    sanitize_times();

    if (seen_time != Mem->LastSightTime)
    {
        // if ( seen_time != Mem->LastSightTime-1 ){
        my_error("Why the sight delay? %d %d (%d %d) type %d %c %d",
                 seen_time.t, Mem->LastSightTime.t, seen, seen_moving, type,
                 ((PlayerObject *)this)->side, ((PlayerObject *)this)->unum);
        //}
        return;
    }

    if (Mem->MyUpdateTime() != time)
        my_error("Must have missed a cycle (mobile)");
    if (Mem->LastSightTime != Mem->CurrentTime && Mem->LastSightTime != Mem->CurrentTime - 1)
        my_error("Didn't see in the past cycle");

    /****** Position ********/

    if (seen)
    {
        /* Update gpos from dist,dir */
        if (patime != time || pdtime != time)
        {
            my_error("Need to do something here--got just ang %d %d  %d %d",
                     pdtime.t, pdtime.s, time.t, time.s); /*skip update if no dist */
        }
        rnpos = Polar2Vector(dist, ang_from_neck);
        rntime = time;
        gpos = Mem->MyPos() + rnpos.rotate(Mem->MyNeckGlobalAng());
        gtime = time;
        gconf = max_conf;
    }

    /****** Velocity ********/

    if (seen_moving)
    {
        /* Update velocity from distch,dirch */
        /* This isn't the same rvel returned by get_rel_vel(), so don't store it */
        /* this way of computing the realtive velocity didn't work all the time
           especially when the ball was approaching you
        Vector temp_rvel = Vector(distch, dist*Tan(dirch)).rotate(ang_from_neck);
        */
        /* printf("rvel:  (%f %f)\n",rvel.mod(),rvel.dir()); */

        /* Now compute the relative velocity the right way and see if they agree */
        Vector temp2_rvel;
        Vector rpos = Polar2Vector(dist, ang_from_neck);
        rpos = rpos.Normalize();
        temp2_rvel.x = distch * rpos.x - (dirch * M_PI / 180 * dist * rpos.y);
        temp2_rvel.y = distch * rpos.y + (dirch * M_PI / 180 * dist * rpos.x);
#ifdef NOT_NEEDED_NOW
        /* now check to see if we get different answers */
        if (fabs(temp_rvel.x - temp2_rvel.x) > .01 ||
            fabs(temp_rvel.y - temp2_rvel.y) > .01)
        {
            // my_error("rpos: (%6.2f, %6.2f)",rpos.x, rpos.y);
            Vector temp_gvel = Mem->MyVel() + temp_rvel.rotate(Mem->MyNeckGlobalAng());
            Vector temp2_gvel = Mem->MyVel() + temp2_rvel.rotate(Mem->MyNeckGlobalAng());
            my_error("Diff vels! old: (%6.2f, %6.2f)   new: (%6.2f, %6.2f)",
                     temp_gvel.x, temp_gvel.y, temp2_gvel.x, temp2_gvel.y);
        }
#endif

        // gvel = Mem->MyVel() + temp_rvel.rotate(Mem->MyNeckGlobalAng());
        gvel = Mem->MyVel() + temp2_rvel.rotate(Mem->MyNeckGlobalAng());
        gvtime = time;
        gvconf = max_conf;
    }
}

/********************************************************************************/

void MobileObject::update_heard(Time time)
{
    if (heard)
    {
        if (time < hptime)
            my_error("How did I fall behind?");
        if (pos_valid() && gtime < time)
            my_error("Should estimate before processing sounds (m1)");

        float decayed_hpconf = hpconf * Exp(conf_decay, (time - hptime));
        if (decayed_hpconf > gconf || !pos_valid() ||
            (decayed_hpconf == gconf && get_dist() > hpdist))
        {
            /* Update gpos from hpos     */
            gpos = hpos;
            gtime = hptime;
            gconf = hpconf;
            if (type == OBJ_Ball)
                Mem->LogAction2(200, "Updating the ball's position based on heard info");
        }
    }

    if (heard_moving)
    {
        if (time < hvtime)
            my_error("How did I fall behind?");
        if (vel_valid() && gvtime < time)
            my_error("Should estimate before processing sounds (m2) %d %d %d %d", gvtime.t, gvtime.s, time.t, time.s);

        float decayed_hvconf = hvconf * Exp(conf_decay, (time - hptime));
        if ((!pos_valid() || get_dist() > Mem->SP_feel_distance || !vel_valid()) &&
            /* Don't listen to vel if object's near */
            (decayed_hvconf > gvconf || !pos_valid() ||
             (decayed_hvconf == gvconf && get_dist() > hvdist)))
        {
            /* Update gvel from hvel */
            gvel = hvel;
            gvtime = hvtime;
            gvconf = hvconf;
            if (type == OBJ_Ball)
                Mem->LogAction2(200, "Updating the ball's velocity based on heard info");
        }
    }

    /* keep updating hptime, hvtime until at time */
    while (vel_valid() && gvtime < gtime)
    {
        estimate_vel(gvtime + 1);
        ++gvtime;
        gvconf *= conf_decay;
    }

    while (pos_valid() && gtime < time)
        update_estimate(gtime + 1);
}

/********************************************************************************/

void MobileObject::update_estimate(Time time)
{
    if (!pos_valid())
        return;

    if (gtime == time && vel_valid())
    { /* just vel */
        if (gvtime == time)
            my_error("pos and vel already updated");
        if (Mem->NewAction &&
            Mem->LastActionValid(gvtime) &&
            Mem->LastActionType() == CMD_kick &&
            type == OBJ_Ball)
            update_kick(gvtime);

        while (vel_valid() && gvtime < gtime)
        {
            estimate_vel(time);
            ++gvtime;
            gvconf *= conf_decay;
        }
        return;
    }

    while (gtime < time)
    { /* both pos and vel */

        if (Mem->NewAction && Mem->LastActionValid(gtime) &&
            Mem->LastActionType() == CMD_kick && type == OBJ_Ball && vel_valid())
            update_kick(gtime);

        estimate_pos(time);
        ++gtime;
        gconf *= conf_decay;

        while (vel_valid() && gvtime < gtime)
        {
            estimate_vel(time);
            ++gvtime;
            gvconf *= conf_decay;
        }
    }
}

/********************************************************************************/

Vector MobileObject::estimate_future_pos(int steps, Vector extra_vel, Vector extra_vel_per)
{
    if (!pos_valid())
        my_error("Can't estimate future if don't know present");

    Vector position = gpos;
    Vector velocity;
    if (vel_valid())
        velocity = gvel + extra_vel;
    else
        velocity = 0;

    for (int i = 0; i < steps; i++)
    {
        velocity += extra_vel_per;
        if (velocity.mod() > max_speed)
            velocity *= (max_speed / velocity.mod());
        position += velocity;
        velocity *= motion_decay;
    }

    return position;
}

/********************************************************************************/

void MobileObject::reset()
{
    Object::reset();
    gconf = 0;
    hpconf = hvconf = gvconf = 0;
    gvtime = hptime = hvtime = 0;
    hpdist = hvdist = 0;
    heard = seen_moving = heard_moving = FALSE;
}

/********************************************************************************/

void MobileObject::clear_seen()
{
    Object::clear_seen();
    seen_moving = FALSE;
}

/********************************************************************************/

void MobileObject::forget()
{
    gconf = gvconf = 0;
#ifndef NO_ACTION_LOG
    if (type == OBJ_Ball)
    {
        Vector pos = get_rel_to_neck_pos();
        Mem->LogAction6(175, "Forgetting ball: (%.1f, %.1f) (%.1f, %.1f)",
                        pos.x, pos.y, (pos - Vector(3, 0)).x, (pos - Vector(3, 0)).y);
        Mem->LogAction4(175, "Forgetting ball (still): old va: %.0f  new va: %.0f",
                        Mem->MyViewAngle(Mem->LastSightTime),
                        Mem->MyViewAngle(Mem->CurrentTime));
    }
#endif
}

/********************************************************************************/

void MobileObject::sanitize_times()
{
    Object::sanitize_times();
    Mem->sanitize_time(hptime);
    Mem->sanitize_time(hvtime);
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        BallObject Class                                      */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void BallObject::Initialize(float max, float min_valid, float decay, float motion, float max_speed)
{
    MobileObject::Initialize(OBJ_Ball, max, min_valid, decay, motion, max_speed);
    reset();
    ektime = 0;
    use_pos_based_vel_estimate = TRUE;
    pos_based_vel_time = 0;
    past_kick_idx = 0;
    for (int i = 0; i < num_past_kicks; i++)
        past_kicks[i].time = -1;
}

/********************************************************************************/

Bool BallObject::moving()
{
    /* experimentally checked that if just the player's moving, ball speed could register > .15 */
    return vel_valid() && get_speed() > Mem->CP_ball_moving_threshold ? TRUE : FALSE;
}

/********************************************************************************/

Bool BallObject::kickable(float buffer)
{
    if (!pos_valid() || !Mem->MyConf())
        return FALSE;
    return get_dist() <= Mem->SP_kickable_area - buffer ? TRUE : FALSE;
}

/********************************************************************************/

Bool BallObject::catchable()
{
    if (!pos_valid() || !Mem->MyConf())
        return FALSE;
    return get_dist() <= Mem->SP_catch_area_l ? TRUE : FALSE;
}

/********************************************************************************/

float BallObject::get_kick_rate(Time time)
{
    if (ektime != time)
    {
        if (!kickable())
        {
            my_error("not kickable %d %d  %d %d", ektime.t, ektime.s, time.t, time.s);
            my_error("not kickable (ball at %.1f %.1f)", gpos.x, gpos.y);
        }

        if (time != Mem->MyUpdateTime())
            my_error("Calculating a kick_rate for a time not equal to update time (%d.%d) (%d.%d)",
                     time.t, time.s, Mem->MyUpdateTime().t, Mem->MyUpdateTime().s);

        effective_kick_rate = calc_kick_rate();
        ektime = time;
    }
    return effective_kick_rate;
}

/********************************************************************************/

void BallObject::update_kick(Time time)
{
    if (!(Mem->LastActionValid(time)))
        my_error("No action at that time (kick)");
    if (Mem->LastActionType() != CMD_kick)
        my_error("Last action wasn't a kick");
    if (!pos_valid() || !vel_valid())
        my_error("Can't update ball if not known before");
    if (!Mem->MyConf())
        my_error("Can't update ball if not localized");

    Mem->LogAction7(180, "Updating kick: time %d.%d, kick %.2f at %.2f, rate %.5f",
                    time.t, time.s, Mem->LastActionPower(), Mem->LastActionAngle(),
                    get_kick_rate(time));

    Vector kick_effect =
        Polar2Vector(get_kick_rate(time) * Mem->LastActionPower(),
                     Mem->MyBodyAng() + Mem->LastActionAngle());
    /* store this kick for later estimation */
    /* this is now done explicitly by send_action in client.C
    past_kicks[past_kick_idx].time = time;
    past_kicks[past_kick_idx].kick_effect = kick_effect;
    past_kick_idx = past_kick_inc(past_kick_idx);
    */

    gvel += kick_effect;
    if (gvel.mod() > max_speed)
        gvel *= (max_speed / gvel.mod());
}

/********************************************************************************/

float BallObject::calc_kick_rate(float dist, float ang)
{
    return Mem->SP_kick_power_rate *
           (1 - .25 * fabs(ang) / 180.0 -
            .25 * (dist - Mem->SP_ball_size - Mem->SP_player_size) / Mem->SP_kickable_margin);
}

/********************************************************************************/

void BallObject::set_past_kick(float pow, AngleDeg ang, Time t)
{
    Vector kick_effect =
        Polar2Vector(get_kick_rate(t) * pow, Mem->MyBodyAng() + ang);
    past_kicks[past_kick_idx].time = t;
    past_kicks[past_kick_idx].kick_effect = kick_effect;
    past_kick_idx = past_kick_inc(past_kick_idx);
}

/********************************************************************************/

void BallObject::forget_past_kick(Time t)
{
    for (int i = 0; i < num_past_kicks; i++)
    {
        if (past_kicks[i].time == t)
            past_kicks[i].kick_effect = Vector(0, 0); /* make the kick effectless */
    }
}

/********************************************************************************/

void BallObject::update(Time time)
{

    MobileObject::update(time);

    if (kickable())
    {
        effective_kick_rate = calc_kick_rate();
        ektime = time;
    }

    if (vel_valid() && get_speed() > max_speed)
        gvel *= max_speed / get_speed();
}

/********************************************************************************/

void BallObject::update_seen(Time time)
{
    Time prev_seen_time = last_seen_time;
    Vector prev_seen_pos = last_seen_pos;

    Vector epos, evel;
    float estimate_valid = pos_valid();
    if (estimate_valid)
    {
        /* if LastSightTime == CurrentTime-1, current is an estimate for the appropriate time */
        epos = gpos;
        evel = gvel;
        if (Mem->LastSightTime == Mem->CurrentTime)
        {
            if (time != Mem->CurrentTime)
                my_error("BallObject::update seen: times seen strange: %d.%d %d.%d",
                         time.t, time.s, Mem->CurrentTime.t, Mem->CurrentTime.s);
            /* see if there is a kick we need to update for */
            if (Mem->LastActionValid(Mem->CurrentTime - 1) &&
                Mem->LastActionType() == CMD_kick &&
                Mem->LastActionTime() == Mem->CurrentTime - 1)
            {
                /*Mem->LogAction5(180, "ball velocity inval: updating for kick %.2f at %.2f, rate %.5f",
                        Mem->LastActionPower(), Mem->LastActionAngle(), get_kick_rate(Mem->CurrentTime-1)); */
                evel += Polar2Vector(get_kick_rate(Mem->CurrentTime - 1) * Mem->LastActionPower(),
                                     Mem->MyBodyAng() + Mem->LastActionAngle());
            }
            epos += evel; /* estimate is the old position + velocity */
            evel *= Mem->SP_ball_decay;
            estimate_valid = vel_valid();
        }
    }

    MobileObject::update_seen(time);

    /*****************************************************/
    /* THIS IS DISTANCE BASED BALL VELOCITY INVALIDATION */
    /*****************************************************/
    /* Only if I see the object, but not its velocity and I have a valid estimate */
    if (seen && !seen_moving && estimate_valid && Mem->sight_position_correction_time == time)
    {
        /* first we need to see if we estiamted a collision, in which case we should
           NOT invalidate the velocity, but multiply it by -1.0 instead */
        if ((epos + Mem->sight_position_correction).dist(Mem->MyPos()) <
            Mem->SP_player_size + Mem->SP_ball_size)
        {
            Mem->LogAction3(175, "ball: update seen predicts a collision %.2f",
                            (epos + Mem->sight_position_correction).dist(Mem->MyPos()));
            gvel = evel * -.1;
            gvconf = Mem->CP_min_valid_conf / Mem->CP_ball_conf_decay; // only valid for 2 cycles
            gvtime = time;
            /* now, don't use postion based velocity estimation */
            use_pos_based_vel_estimate = FALSE;
            pos_based_vel_time = time;
        }
        else
        {
            /* first we need to figure out the maximum error that the server would give
             total_dist_err is supposed to the the total dist across the error margin of
             the reporting of the ball's position. The angle is +/- .5 degrees and the
             distance is quantize as follows
             Quantize(exp(Quantize(log(vi.distance + EPS),vi.qstep)),0.1) ;

             After that, we have to add in the possible noise from the server.
             We're just going to add in the noise once for every cycle between our
             sight times
          */
            float d = get_dist();
            float just_dist_err = d * Mem->quantize_err_const / 2;
            float inner_perp_err = Mem->Tan_of_half_deg * (d - just_dist_err);
            float outer_perp_err = Mem->Tan_of_half_deg * (d + just_dist_err);
            float total_dist_err =
                sqrt(Sqr(inner_perp_err) + Sqr(just_dist_err)) +
                sqrt(Sqr(outer_perp_err) + Sqr(just_dist_err));

            total_dist_err += Mem->LastSightInterval *
                              sqrt(2.0 * Sqr(Mem->SP_ball_rand * gvel.mod()));

            Vector diff = gpos - (epos + Mem->sight_position_correction);
            if (diff.mod() > total_dist_err * Mem->CP_ball_vel_invalidation_factor)
            {
                Mem->LogAction6(175, "Invalidating ball velocity: %.2f > %.2f, thought vel was (%.2f, %.2f)",
                                diff.mod(), total_dist_err * Mem->CP_ball_vel_invalidation_factor,
                                evel.x, evel.y);
                Mem->LogAction6(175, "Invalidating ball vel(still): gpos (%.1f %.1f), sight_position_correction (%.1f %.1f)",
                                gpos.x, gpos.y,
                                Mem->sight_position_correction.x, Mem->sight_position_correction.y);
                gvconf = 0;
                if (!(pos_based_vel_time == time && !use_pos_based_vel_estimate))
                {
                    // if it's already been set, don't touch it
                    use_pos_based_vel_estimate = (time - prev_seen_time == 1) ? TRUE : FALSE;
                    pos_based_vel_time = time;
                }
            }
        }
    }

    /**********************************************/
    /* THIS IS POSITION BASED VELOCITY ESTIMATION */
    /**********************************************/
    if (Mem->MyConf() && pos_valid() &&
        gvconf <= Exp(Mem->CP_ball_conf_decay, 2) &&
        get_dist() <= Mem->SP_feel_distance &&
        prev_seen_time >= Mem->PreviousSightTime() &&
        !(use_pos_based_vel_estimate == FALSE &&
          pos_based_vel_time == time))
    {
        /* Don't estimate velocity if the ball's far -- too much noise */

        /* Can estimate based on the last seen position of the ball */
        // cout << "Time: " << time.t << "  Using position based vel estimate" << endl;
        Vector total_kick_eff(0, 0);
        /* we're goign to look through our past kicks to find ones that occured in
           the time cycle that we're looking at */
        int pkidx;
        for (pkidx = past_kick_dec(past_kick_idx);
             past_kicks[pkidx].time >= prev_seen_time;
             pkidx = past_kick_dec(pkidx))
        {
            /* this kick falls in our range */
            if (past_kicks[pkidx].time > time)
                my_error("Postion Based Vel Estimate: Already have a future kick???");
            total_kick_eff += past_kicks[pkidx].kick_effect *
                              SumGeomSeries(1, motion_decay, time - past_kicks[pkidx].time);
        }
        /* POSITION BASED VELOCITY ESTIMATION
           A problem that arises using position based velocity estimes is that
           we store the global position of the ball as calculated relative to the
           player (whose global position we identify by flags). Becuase of the
           error in seeing the flags, our global position ossilates around the true
           value, adding additional error to the position based velocity
           calculation. The solution is to take the difference between where we
           expected to be and where we observe we are at each new sight and store
           that. We then correct the global ball position for that. This
           essentially removes the error in calculating our global position from
           this calculation, leaving only the ball observation error, giving us a
           better value */
        if (Mem->CP_use_new_position_based_vel && Mem->MyConf() &&
            Mem->sight_position_correction_time == time)
        {
            // cout << "time: " << time.t << "\tcorr: " << Mem->sight_position_correction.mod() << endl;
            gvel = (gpos - (prev_seen_pos + Mem->sight_position_correction) - total_kick_eff) / SumGeomSeries(1, motion_decay, time - prev_seen_time);
            Mem->LogAction8(175, "Position based velocity estimating: gpos (%.1f %.1f), prev_seen_pos (%.1f %.1f), sight_position_correction (%.1f %.1f)",
                            gpos.x, gpos.y, prev_seen_pos.x, prev_seen_pos.y,
                            Mem->sight_position_correction.x, Mem->sight_position_correction.y);
        }
        else
        {
            Mem->LogAction6(175, "Old position based velocity estimating: gpos (%.1f %.1f), prev_seen_pos (%.1f %.1f)",
                            gpos.x, gpos.y, prev_seen_pos.x, prev_seen_pos.y);
            gvel = (gpos - prev_seen_pos - total_kick_eff) / SumGeomSeries(1, motion_decay, time - prev_seen_time);
            /* gvel is now the velocty at prev_seen_time */
        }

        Time t;
        pkidx = past_kick_inc(pkidx);
        for (t = prev_seen_time; t < time; ++t)
        {
            if (pkidx != past_kick_idx && past_kicks[pkidx].time == t)
            {
                gvel += past_kicks[pkidx].kick_effect;
                pkidx = past_kick_inc(pkidx);
            }
            gvel *= motion_decay;
        }

        gvconf = max_conf; /* so it overrides heard info */
        gvtime = time;
    }

    if (seen_time > last_seen_time)
    {
        last_seen_time = time; /* == last_seen_time except when clock was stopped */
        last_seen_pos = gpos;
    }
}

/********************************************************************************/

void BallObject::estimate_pos(Time time)
{
    MobileObject::estimate_pos(time);
    if (!Mem->MyConf())
        return;
    if (!pos_valid())
        return; /* Can't check for collisions */

    /* Only worry about collisions for the ball */
    if (get_dist() < Mem->SP_player_size + Mem->SP_ball_size)
    {
        if (vel_valid())
        {
            float r = Mem->SP_ball_size + Mem->SP_player_size;
            float d = get_dist();
            float th = fabs(GetNormalizeAngleDeg(get_rel_to_body_heading() - 180));
            float l1 = d * Cos(th);
            float h = d * Sin(th);
            float cosp = h / r;
            float sinp = sqrt(1.0 - Sqr(cosp));
            float l2 = r * sinp;
            gpos += gvel * motion_decay * (-(l1 + l2) / Max(get_speed() * motion_decay, 1.0e-10));
            gvel *= -0.1;
            Mem->LogAction2(160, "Ball collision");
            /* turn off position based velocity estimation for this cycle */
            use_pos_based_vel_estimate = FALSE;
            pos_based_vel_time = time;
        }
        /* my_stamp; printf("COLLISION!  --   check computation\n"); */
    }
}

/********************************************************************************/

Bool BallObject::in_view_range(AngleDeg view_ang, float angle_buffer, float distance_buffer)
{
    if (!pos_valid() || !Mem->MyConf())
        return FALSE;

    if (MobileObject::in_view_range(view_ang, angle_buffer, 0) ||
        get_dist() < Mem->SP_feel_distance - distance_buffer)
        return TRUE;
    return FALSE;
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        PlayerObject Class                                    */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void PlayerObject::Initialize(float max, float min_valid, float decay, float motion, float max_speed)
{
    MobileObject::Initialize(OBJ_Player, max, min_valid, decay, motion, max_speed);
    side = 'f'; /* free */
    unum = Unum_Unknown;
    reset();
}

/********************************************************************************/

AngleDeg PlayerObject::get_rel_to_body_body_ang()
{
    float f = get_abs_body_ang() - Mem->MyBodyAng();
    NormalizeAngleDeg(&f);
    return f;
}

/********************************************************************************/

AngleDeg PlayerObject::get_rel_to_body_neck_ang()
{
    float f = get_abs_neck_ang() - Mem->MyBodyAng();
    NormalizeAngleDeg(&f);
    return f;
}

/********************************************************************************/

AngleDeg PlayerObject::get_rel_to_neck_body_ang()
{
    float f = get_abs_body_ang() - Mem->MyNeckGlobalAng();
    NormalizeAngleDeg(&f);
    return f;
}

/********************************************************************************/

AngleDeg PlayerObject::get_rel_to_neck_neck_ang()
{
    float f = get_abs_neck_ang() - Mem->MyNeckGlobalAng();
    NormalizeAngleDeg(&f);
    return f;
}

/********************************************************************************/

AngleDeg PlayerObject::get_neck_rel_ang()
{
    float f = get_rel_to_neck_neck_ang() - get_rel_to_neck_body_ang();
    NormalizeAngleDeg(&f);
    return f;
}

/********************************************************************************/

void PlayerObject::set_body_ang_from_neck(AngleDeg body_ang, Time time)
{
    rnbang = body_ang;
    rnbatime = time;
    seen_body_ang = TRUE;
}

/********************************************************************************/

void PlayerObject::set_neck_ang_from_neck(AngleDeg neck_ang, Time time)
{
    rnnang = neck_ang;
    rnnatime = time;
    seen_neck_ang = TRUE;
}

/********************************************************************************/

void PlayerObject::set_heard_info_w_angs(float x, float y, float pconf, float dx, float dy, float vconf,
                                         AngleDeg bang, float bconf, AngleDeg nang, float nconf,
                                         float dist, Time time)
{
    /* Don't overwrite better infor from the same cycle      */
    /* If I have better face info, I have better pos. info   */
    /* If I have better face info, I have better neck info   */
    /* IMPORTANT -- above assumes body and neck angles are only heard together
       and always with position info */
    if (heard_body_ang &&
        (bconf < hbaconf || (bconf == hbaconf && dist > hbadist)))
        return;

    MobileObject::set_heard_info(x, y, pconf, dx, dy, vconf, dist, time);
    hbang = bang;
    hbaconf = bconf;
    hbadist = dist;
    hbatime = time;
    heard_body_ang = TRUE;

    hnang = nang;
    hnaconf = nconf;
    hnadist = dist;
    hnatime = time;
    heard_neck_ang = TRUE;
}

/********************************************************************************/

void PlayerObject::update(Time time)
{
    /* IMPORTANT -- assuming that body and neck angle are only seen together
       Otherwise, might have to reason about whether the neck angle has changed */

    if (Mem->NewSight && seen_body_ang)
    {
        if (!seen_neck_ang)
            my_error("Should see body and neck angle together");
        update_seen_body_and_neck_angs(Mem->LastSightTime);
    }

    if (gbatime < time)
    {
        if (gnatime >= time)
            my_error("Should see body and neck angle together");
        update_estimate(time);
    }

    if (heard_body_ang)
    {
        if (!heard_neck_ang)
            my_error("Should hear body and neck angle together");
        update_heard(time);
    }

    MobileObject::update(time);
    seen_body_ang = heard_body_ang = seen_neck_ang = heard_neck_ang = FALSE;
}

/********************************************************************************/

void PlayerObject::update_seen_body_and_neck_angs(Time time)
{
    if (!Mem->MyConf())
        return;

    sanitize_times();

    if (seen_time != Mem->LastSightTime)
    {
        /* if ( seen_time != Mem->LastSightTime-1 ) */
        my_error("Why the sight delay(2)? %d (%d %d) %d %d %c %d",
                 seen_time.t, seen_body_ang, seen_neck_ang,
                 rnbatime.t, rnnatime.t, side, unum);
        return;
    }

    if (Mem->MyUpdateTime() != time)
        my_error("Must have missed a cycle (player)");

    gbang = GetNormalizeAngleDeg(rnbang + Mem->MyNeckGlobalAng());
    gbaconf = max_conf;
    gbatime = time;
    gnang = GetNormalizeAngleDeg(rnnang + Mem->MyNeckGlobalAng());
    /**** OR gbang + rnnang ****/
    gnaconf = max_conf;
    gnatime = time;
}

/********************************************************************************/

void PlayerObject::update_estimate(Time time)
{
    if (!body_ang_valid() || !neck_ang_valid())
        return;

    while (gbatime < time)
    {
        gbaconf *= conf_decay;
        ++gbatime;
    }
    while (gnatime < time)
    {
        gnaconf *= conf_decay;
        ++gnatime;
    }
}

/********************************************************************************/

void PlayerObject::update_heard(Time time)
{
    if (!heard_body_ang || !heard_neck_ang) /* Should only be together */
        my_error("Should only be here if I heard face info");
    if (time < hbatime || time < hnatime)
        my_error("How did I fall behind?");
    if (body_ang_valid() && gbatime < time)
        my_error("Should estimate before processing sounds (p1)");
    if (neck_ang_valid() && gnatime < time)
        my_error("Should estimate before processing sounds (p2)");

    float decayed_hbaconf = hbaconf * Exp(conf_decay, (time - hbatime));
    if (decayed_hbaconf > gbaconf || !pos_valid() ||
        (decayed_hbaconf == gbaconf && get_dist() > hbadist))
    {
        /* Update gpos from hpos     */
        gbang = hbang;
        gbatime = hbatime;
        gbaconf = hbaconf;
    }

    while (body_ang_valid() && gbatime < time)
        update_estimate(gbatime + 1);

    float decayed_hnaconf = hnaconf * Exp(conf_decay, (time - hnatime));
    if (decayed_hnaconf > gnaconf || !pos_valid() ||
        (decayed_hnaconf == gnaconf && get_dist() > hnadist))
    {
        /* Update gpos from hpos     */
        gnang = hnang;
        gnatime = hnatime;
        gnaconf = hnaconf;
    }

    while (neck_ang_valid() && gnatime < time)
        update_estimate(gnatime + 1);
}

/********************************************************************************/

void PlayerObject::reset()
{
    MobileObject::reset();
    hbaconf = gbaconf = hnaconf = gnaconf = 0;
    rnbatime = gbatime = hbatime = rnnatime = gnatime = hnatime = 0;
    hbadist = hnadist = 0;
    seen_body_ang = heard_body_ang = seen_neck_ang = heard_neck_ang = FALSE;
}

/********************************************************************************/

void PlayerObject::clear_seen()
{
    MobileObject::clear_seen();
    seen_body_ang = seen_neck_ang = FALSE;
}

/********************************************************************************/

void PlayerObject::forget()
{
    MobileObject::forget();
    gbaconf = gnaconf = 0;
#ifndef NO_ACTION_LOG
    float slide_dist = 5 / Sin(Mem->MyViewAngle(Mem->CurrentTime));
    Vector pos = get_rel_to_neck_pos();
    Mem->LogAction7(175, "Forgetting player %d: (%.1f, %.1f) (%.1f, %.1f)",
                    (side == Mem->MySide) ? unum : -unum,
                    pos.x, pos.y,
                    (pos - Vector(slide_dist, 0)).x, (pos - Vector(slide_dist, 0)).y);
    Mem->LogAction5(175, "Forgetting player %d (still): old va: %.0f  new va: %.0f",
                    (side == Mem->MySide) ? unum : -unum,
                    Mem->MyViewAngle(Mem->LastSightTime),
                    Mem->MyViewAngle(Mem->CurrentTime));
#endif
}

/********************************************************************************/

void PlayerObject::sanitize_times()
{
    MobileObject::sanitize_times();
    Mem->sanitize_time(rnbatime);
    Mem->sanitize_time(hbatime);
    Mem->sanitize_time(rnnatime);
    Mem->sanitize_time(hnatime);
}

/********************************************************************************/

Bool PlayerObject::in_view_range(AngleDeg view_ang, float angle_buffer, float distance_buffer)
{
    if (!pos_valid() || !Mem->MyConf())
        return FALSE;

    if (MobileObject::in_view_range(view_ang, angle_buffer, 0) &&
        get_dist() < Mem->SP_unum_far_length - distance_buffer)
        return TRUE;
    return FALSE;
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
/*                        PositionInfo Class                                    */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void PositionInfo::Initialize()
{
    /* printf("Calling Position Initialize\n"); */

    /* if true, multiply all coords by -1 inside initialize */
    Bool rotate = (Mem->MySide == 'l') ? FALSE : TRUE;

    int i = 0;

    Marker = new StationaryObject[SP_num_markers];
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0, 0.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Goal_L */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0, 0.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Goal_R */

    Marker[i].Initialize((MarkerType)i, Vector(0.0, 0.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_C */
    Marker[i].Initialize((MarkerType)i, Vector(0.0, -SP_pitch_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_CT */
    Marker[i].Initialize((MarkerType)i, Vector(0.0, SP_pitch_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_CB */

    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0, -SP_pitch_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_LT */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0, SP_pitch_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_LB */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0, -SP_pitch_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_RT */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0, SP_pitch_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_RB */

    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 + SP_penalty_area_length, -SP_penalty_area_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_PLT */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 + SP_penalty_area_length, 0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_PLC */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 + SP_penalty_area_length, SP_penalty_area_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_PLB */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 - SP_penalty_area_length, -SP_penalty_area_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_PRT */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 - SP_penalty_area_length, 0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_PRC */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 - SP_penalty_area_length, SP_penalty_area_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_PRB */

    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0, -SP_goal_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_GLT */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0, SP_goal_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_GLB */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0, -SP_goal_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_GRT */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0, SP_goal_width / 2.0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_GRB */

    Marker[i].Initialize((MarkerType)i, Vector(-50.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TL50 */
    Marker[i].Initialize((MarkerType)i, Vector(-40.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TL40 */
    Marker[i].Initialize((MarkerType)i, Vector(-30.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TL30 */
    Marker[i].Initialize((MarkerType)i, Vector(-20.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TL20 */
    Marker[i].Initialize((MarkerType)i, Vector(-10.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TL10 */
    Marker[i].Initialize((MarkerType)i, Vector(0.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_T0 */
    Marker[i].Initialize((MarkerType)i, Vector(10.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TR10 */
    Marker[i].Initialize((MarkerType)i, Vector(20.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TR20 */
    Marker[i].Initialize((MarkerType)i, Vector(30.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TR30 */
    Marker[i].Initialize((MarkerType)i, Vector(40.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TR40 */
    Marker[i].Initialize((MarkerType)i, Vector(50.0, -SP_pitch_width / 2.0 - SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_TR50 */

    Marker[i].Initialize((MarkerType)i, Vector(-50.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BL50 */
    Marker[i].Initialize((MarkerType)i, Vector(-40.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BL40 */
    Marker[i].Initialize((MarkerType)i, Vector(-30.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BL30 */
    Marker[i].Initialize((MarkerType)i, Vector(-20.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BL20 */
    Marker[i].Initialize((MarkerType)i, Vector(-10.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BL10 */
    Marker[i].Initialize((MarkerType)i, Vector(0.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_B0 */
    Marker[i].Initialize((MarkerType)i, Vector(10.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BR10 */
    Marker[i].Initialize((MarkerType)i, Vector(20.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BR20 */
    Marker[i].Initialize((MarkerType)i, Vector(30.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BR30 */
    Marker[i].Initialize((MarkerType)i, Vector(40.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BR40 */
    Marker[i].Initialize((MarkerType)i, Vector(50.0, SP_pitch_width / 2.0 + SP_pitch_margin),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_BR50 */

    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 - SP_pitch_margin, -30),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_LT30 */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 - SP_pitch_margin, -20),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_LT20 */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 - SP_pitch_margin, -10),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_LT10 */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 - SP_pitch_margin, 0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_L0 */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 - SP_pitch_margin, 10),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_LB10 */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 - SP_pitch_margin, 20),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_LB20 */
    Marker[i].Initialize((MarkerType)i, Vector(-SP_pitch_length / 2.0 - SP_pitch_margin, 30),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_LB30 */

    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 + SP_pitch_margin, -30),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_RT30 */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 + SP_pitch_margin, -20),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_RT20 */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 + SP_pitch_margin, -10),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_RT10 */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 + SP_pitch_margin, 0),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_R0 */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 + SP_pitch_margin, 10),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_RB10 */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 + SP_pitch_margin, 20),
                         CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* Flag_RB20 */
    Marker[i].Initialize((MarkerType)i, Vector(SP_pitch_length / 2.0 + SP_pitch_margin, 30),
                         CP_max_conf, CP_min_valid_conf, rotate); /* Flag_RB30 */

    if (MySide == 'l')
    {
        RM_My_Goal = Goal_L;
        RM_Their_Goal = Goal_R;
        RM_LB_Flag = Flag_LT;
        RM_LC_Flag = Flag_CT;
        RM_LF_Flag = Flag_RT;
        RM_RB_Flag = Flag_LB;
        RM_RC_Flag = Flag_CB;
        RM_RF_Flag = Flag_RB;
        RM_My_PC_Flag = Flag_PLC;    /* Center of my penalty area */
        RM_Their_PC_Flag = Flag_PRC; /* Center of theirs          */
    }
    else
    { /* MySide == 'r' */
        RM_My_Goal = Goal_R;
        RM_Their_Goal = Goal_L;
        RM_LB_Flag = Flag_RB;
        RM_LC_Flag = Flag_CB;
        RM_LF_Flag = Flag_LB;
        RM_RB_Flag = Flag_RT;
        RM_RC_Flag = Flag_CT;
        RM_RF_Flag = Flag_LT;
        RM_My_PC_Flag = Flag_PRC;    /* Center of my penalty area */
        RM_Their_PC_Flag = Flag_PLC; /* Center of theirs          */
    }

    i = 0;

    Fieldline = new StationaryObject[SP_num_lines];
    Fieldline[i].Initialize((SideLine)i, Vector(-SP_pitch_length / 2.0, 0.0),
                            CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* SL_Left */
    Fieldline[i].Initialize((SideLine)i, Vector(SP_pitch_length / 2.0, 0.0),
                            CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* SL_Right */
    Fieldline[i].Initialize((SideLine)i, Vector(0.0, -SP_pitch_width / 2.0),
                            CP_max_conf, CP_min_valid_conf, rotate);
    i++; /* SL_Top */
    Fieldline[i].Initialize((SideLine)i, Vector(0.0, SP_pitch_width / 2.0),
                            CP_max_conf, CP_min_valid_conf, rotate); /* SL_Bottom */

    Ball.Initialize(CP_max_conf, CP_min_valid_conf, CP_ball_conf_decay, SP_ball_decay, SP_ball_speed_max);

    num_players = SP_team_size * 2 - 1;

    Player = new PlayerObject[num_players]; /* allow for all players but me          */
    for (i = 0; i < num_players; i++)
        Player[i].Initialize(CP_max_conf, CP_min_valid_conf, CP_player_conf_decay,
                             SP_player_decay, SP_player_speed_max);

    for (i = 0; i < num_players; i++)
        FreePlayer[i] = &(Player[(num_players - 1) - i]); /* Player array backwards: take from end */

    UnknownPlayer = new TempPlayerObject[num_players];

    num_seen_markers = 0;
    num_my_players = 0;
    num_their_players = 0;
    num_teamless_players = 0;
    num_free_players = num_players;
    num_unknown_players = 0;

    ClosestMarker = ClosestMotionMarker = No_Marker;
    SeenLine = SL_No_Line;

    for (i = 1; i <= SP_team_size; i++)
        TiredTimes[i] = -CP_say_tired_interval;

    OwnPenaltyArea = Rectangle(-SP_pitch_length / 2,                          // left
                               -SP_pitch_length / 2 + SP_penalty_area_length, // right
                               -SP_penalty_area_width / 2,                    // top
                               SP_penalty_area_width / 2);                    // bottom
    OwnGoalieArea = Rectangle(-SP_pitch_length / 2,                           // left
                              -SP_pitch_length / 2 + SP_goal_area_length,     // right
                              -SP_goal_area_width / 2,                        // top
                              SP_goal_area_width / 2);                        // bottom

    TheirPenaltyArea = Rectangle(SP_pitch_length / 2 - SP_penalty_area_length, // left
                                 SP_pitch_length / 2,                          // right
                                 -SP_penalty_area_width / 2,                   // top
                                 SP_penalty_area_width / 2);                   // bottom
    TheirGoalieArea = Rectangle(SP_pitch_length / 2 - SP_goal_area_length,     // left
                                SP_pitch_length / 2,                           // right
                                -SP_goal_area_width / 2,                       // top
                                SP_goal_area_width / 2);                       // bottom

    FieldRectangle = Rectangle(-SP_pitch_length / 2,
                               SP_pitch_length / 2,
                               -SP_pitch_width / 2,
                               SP_pitch_width / 2);

    MyLeftGoalKickSpot = OwnGoalieArea.TopRightCorner();
    MyRightGoalKickSpot = OwnGoalieArea.BottomRightCorner();
    TheirLeftGoalKickSpot = TheirGoalieArea.TopLeftCorner();
    TheirRightGoalKickSpot = TheirGoalieArea.BottomLeftCorner();

    my_offside_line = 0;

    quantize_err_const = exp(SP_dist_qstep / 2) - exp(-SP_dist_qstep / 2);
    Tan_of_half_deg = Tan(1.0 / 2.0);
}

PositionInfo::~PositionInfo()
{
    delete[] Fieldline;
    delete[] Marker;
    delete[] Player;
    delete[] UnknownPlayer;
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void PositionInfo::SeeLine(SideLine l, float dist, float ang, Time tm)
{
    Fieldline[l].set_polar_from_neck(dist, ang, tm);
    SeenLine = l;
}

/********************************************************************************/

void PositionInfo::SeeLine(SideLine l, float ang, Time tm)
{
    Fieldline[l].set_angle_from_neck(ang, tm);
    SeenLine = l;
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void PositionInfo::SeeMarker(MarkerType marker, float dist, float ang, Time tm)
{
    Marker[marker].set_polar_from_neck(dist, ang, tm);
    SeenMarker[num_seen_markers++] = marker;
}

/********************************************************************************/

void PositionInfo::SeeMarker(MarkerType marker, float ang, Time tm)
{
    Marker[marker].set_angle_from_neck(ang, tm);
    SeenMarker[num_seen_markers++] = marker;
    my_error("Shouldn't process markers when using low quality -- no info");
}

/********************************************************************************/

void PositionInfo::SeeMarker(MarkerType marker, float dist, float ang,
                             float distChng, float dirChng, Time tm)
{
    Marker[marker].set_chinfo(distChng, dirChng, tm);
    Marker[marker].set_polar_from_neck(dist, ang, tm);
    SeenMarker[num_seen_markers++] = marker;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

void PositionInfo::SeeBall(float ang, Time tm)
{
    Ball.set_angle_from_neck(ang, tm);
}

/*********************************************************************************/

void PositionInfo::SeeBall(float dist, float ang, Time tm)
{
    Ball.set_polar_from_neck(dist, ang, tm);
}

/*********************************************************************************/

void PositionInfo::SeeBall(float dist, float ang, float distChng, float dirChng, Time tm)
{
    Ball.set_chinfo(distChng, dirChng, tm);
    SeeBall(dist, ang, tm);
}

/*********************************************************************************/
/********************************************************************************/
/*********************************************************************************/

void PositionInfo::SeePlayer(char side, Unum num, float ang, Time time)
{
    PlayerObject *player;
    if ((player = GetPlayer(side, num)) == NULL)
        player = GetNewPlayer(side, num);

    if (player != NULL)
        player->set_angle_from_neck(ang, time);
    else
        my_error("Can't get a player to see (1) %c %d (%d teamless)", side, num, num_teamless_players);
}

/********************************************************************************/

void PositionInfo::SeePlayer(char side, Unum num, float dist, float ang, Time time)
{
    PlayerObject *player;
    if ((player = GetPlayer(side, num)) == NULL)
        player = GetNewPlayer(side, num);

    if (player != NULL)
        player->set_polar_from_neck(dist, ang, time);
    else
        my_error("Can't get a player to see (2) %c %d (%d teamless)", side, num, num_teamless_players);
}

/********************************************************************************/

void PositionInfo::SeePlayer(char side, Unum num, float dist, float ang,
                             float distChng, float dirChng, float bodydir, float neckdir, Time time)
{
    PlayerObject *player;
    if ((player = GetPlayer(side, num)) == NULL)
        player = GetNewPlayer(side, num);

    if (player != NULL)
    {
        player->set_chinfo(distChng, dirChng, time);
        player->set_body_and_neck_ang_from_neck(bodydir, neckdir, time);
        player->set_polar_from_neck(dist, ang, time);
    }
    else
        my_error("Can't get a player to see (3) %c %d (%d teamless)", side, num, num_teamless_players);
}

/********************************************************************************/

void PositionInfo::SeePlayer(char side, float dist, float ang, Time time)
{
    if (num_unknown_players < num_players)
    {
        UnknownPlayer[num_unknown_players].set(side, dist, ang, time);
        num_unknown_players++;
    }
    else
        my_error("Too many unknown players (1) %d %d %d %d", Mem->LastSightTime.t, Mem->LastSightTime.s,
                 time.t, time.s);
}

/********************************************************************************/

void PositionInfo::SeePlayer(char side, float ang, Time time)
{
    my_error("What can I do that's useful with just the player angle?");
    return;

    if (num_unknown_players < num_players)
    {
        UnknownPlayer[num_unknown_players].set(side, 0, ang, time);
        num_unknown_players++;
    }
    else
        my_error("Too many unknown players (2)");
}

/********************************************************************************/

void PositionInfo::SeePlayer(float dist, float ang, Time time)
{
    if (num_unknown_players < num_players)
    {
        UnknownPlayer[num_unknown_players].set('?', dist, ang, time);
        num_unknown_players++;
    }
    else
        my_error("Too many unknown players (3) %d %d %d %d", Mem->LastSightTime.t, Mem->LastSightTime.s,
                 time.t, time.s);
}

/********************************************************************************/

void PositionInfo::SeePlayer(float ang, Time time)
{
    my_error("What can I do that's useful with just the player angle?");
    return;

    if (num_unknown_players < num_players)
    {
        UnknownPlayer[num_unknown_players].set('?', 0, ang, time);
        num_unknown_players++;
    }
    else
        my_error("Too many unknown players (4)");
}

/********************************************************************************/

void PositionInfo::HearBall(float x, float y, float conf, float dist, Time time)
{
    Ball.set_heard_info(x, y, conf, dist, time);
}

/********************************************************************************/

void PositionInfo::HearBall(float x, float y, float pconf, float dx, float dy, float vconf,
                            float dist, Time time)
{
    Ball.set_heard_info(x, y, pconf, dx, dy, vconf, dist, time);
}

/********************************************************************************/

void PositionInfo::HearPlayer(char side, Unum num, float x, float y, float conf, float dist, Time time)
{
    /* When hearing a player location, remember not to call just getplayer, but also
       closestPlayerto(side,vector) -- might be one of the unknowns
       can just update the side---update_players takes care of the rest */

    PlayerObject *player;
    if ((player = GetPlayer(side, num)) == NULL)
    {
        Vector pos = Vector(x, y);
        if ((player = ClosestPlayerObjectTo(side, pos)) != NULL)
        {
            if (player->side == '?')
                player->side = side;
            player->unum = num;
        }
        else
            player = GetNewPlayer(side, num);
    }

    if (player != NULL)
        player->set_heard_info(x, y, conf, dist, time);
    else
    {
        /* Be more lenient in ClosestPlayerObjectTo??? */
        Unum teammate = ClosestTeammateTo(Vector(x, y));
        Unum opponent = ClosestOpponentTo(Vector(x, y));
        my_error("Can't get a player to hear (1) %c %d (%.1f %.1f)", side, num,
                 teammate == Unum_Unknown ? -1.0 : TeammateDistanceTo(teammate, Vector(x, y)),
                 opponent == Unum_Unknown ? -1.0 : OpponentDistanceTo(opponent, Vector(x, y)));
    }
}

/********************************************************************************/

void PositionInfo::HearPlayer(char side, Unum num, float x, float y, float pconf,
                              float dx, float dy, float vconf, float dist, Time time)
{
    /* When hearing a player location, remember not to call just getplayer, but also
       closestPlayerto(side,vector) -- might be one of the unknowns
       can just update the side---update_players takes care of the rest */

    PlayerObject *player;
    if ((player = GetPlayer(side, num)) == NULL)
    {
        Vector pos = Vector(x, y);
        if ((player = ClosestPlayerObjectTo(side, pos)) != NULL)
        {
            if (player->side == '?')
                player->side = side;
            player->unum = num;
        }
        else
            player = GetNewPlayer(side, num);
    }

    if (player != NULL)
        player->set_heard_info(x, y, pconf, dx, dy, vconf, dist, time);
    else
    {
        /* Be more lenient in ClosestPlayerObjectTo??? */
        Unum teammate = ClosestTeammateTo(Vector(x, y));
        Unum opponent = ClosestOpponentTo(Vector(x, y));
        my_error("Can't get a player to hear (2) %c %d (%.1f %.1f)", side, num,
                 teammate == Unum_Unknown ? -1.0 : TeammateDistanceTo(teammate, Vector(x, y)),
                 opponent == Unum_Unknown ? -1.0 : OpponentDistanceTo(opponent, Vector(x, y)));
    }
}

/********************************************************************************/

void PositionInfo::HearPlayer(char side, Unum num, float x, float y, float pconf, float dx, float dy,
                              float vconf, AngleDeg bang, float bconf, AngleDeg nang, float nconf,
                              float dist, Time time)
{
    /* When hearing a player location, remember not to call just getplayer, but also
       closestPlayerto(side,vector) -- might be one of the unknowns
       can just update the side---update_players takes care of the rest */

    PlayerObject *player;
    if ((player = GetPlayer(side, num)) == NULL)
    {
        Vector pos = Vector(x, y);
        if ((player = ClosestPlayerObjectTo(side, pos)) != NULL)
        {
            if (player->side == '?')
                player->side = side;
            player->unum = num;
        }
        else
            player = GetNewPlayer(side, num);
    }

    if (player != NULL)
        player->set_heard_info_w_angs(x, y, pconf, dx, dy, vconf, bang, bconf, nang, nconf, dist, time);
    else
    {
        /* Be more lenient in ClosestPlayerObjectTo??? */
        Unum teammate = ClosestTeammateTo(Vector(x, y));
        Unum opponent = ClosestOpponentTo(Vector(x, y));
        my_error("Can't get a player to hear (3) %c %d (%.1f %.1f)", side, num,
                 teammate == Unum_Unknown ? -1.0 : TeammateDistanceTo(teammate, Vector(x, y)),
                 opponent == Unum_Unknown ? -1.0 : OpponentDistanceTo(opponent, Vector(x, y)));
    }
}

/********************************************************************************/
/*********************************************************************************/
/********************************************************************************/

float PositionInfo::PlayerPositionValid(char side, Unum n)
{
    if (side == MySide && n == MyNumber)
        return MyConf();

    PlayerObject *player = GetPlayer(side, n);
    if (player == NULL)
        return 0;
    else
        return player->pos_valid();
}

/********************************************************************************/

float PositionInfo::PlayerVelocityValid(char side, Unum n)
{
    if (side == MySide && n == MyNumber)
        return MyVelConf();

    PlayerObject *player = GetPlayer(side, n);
    if (player == NULL)
        return 0;
    else
        return player->vel_valid();
}

/********************************************************************************/

float PositionInfo::PlayerBodyAngleValid(char side, Unum n)
{
    if (side == MySide && n == MyNumber)
        return MyConf();

    PlayerObject *player = GetPlayer(side, n);
    if (player == NULL)
        return 0;
    else
        return player->body_ang_valid();
}

/********************************************************************************/

float PositionInfo::PlayerNeckAngleValid(char side, Unum n)
{
    if (side == MySide && n == MyNumber)
        return MyConf();

    PlayerObject *player = GetPlayer(side, n);
    if (player == NULL)
        return 0;
    else
        return player->neck_ang_valid();
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

Vector PositionInfo::BallPredictedPosition(int steps, float kick_power, AngleDeg kick_angle)
{
    /* kick_angle relative to body */
    if (!BallVelocityValid())
        my_error("Need to know ball velocity for prediction");

    Vector new_vel = Polar2Vector(BallKickRate() * kick_power, MyBodyAng() + kick_angle);
    return GetBall()->estimate_future_pos(steps, new_vel);
}

/********************************************************************************/

Vector PositionInfo::BallPredictedPositionWithQueuedActions(int steps)
{
    if (Action->valid() && Action->type == CMD_kick)
        return BallPredictedPosition(steps, Action->power, Action->angle);
    else
        return BallPredictedPosition(steps);
}

/********************************************************************************/

Bool PositionInfo::BallWillBeKickable(int steps, float dash_power, float buffer)
{
    if (!MyConf() || !BallPositionValid())
        my_error("Can't predict kickable");

    Vector ball_predicted_position = BallPredictedPosition(steps);
    Vector my_predicted_position = MyPredictedPosition(steps, dash_power);

    return my_predicted_position.dist(ball_predicted_position) <= SP_kickable_area - buffer ? TRUE : FALSE;
}

/*****************************************************************************************/

Unum PositionInfo::PlayerWithBall(float buffer)
{
    if (!BallPositionValid())
        return Unum_Unknown;
    Vector ball = BallAbsolutePosition();

    Unum teammate = ClosestTeammateToBall();
    Unum opponent = ClosestOpponentToBall();

    float teammate_distance = 400, opponent_distance = 400;
    if (teammate != Unum_Unknown)
        teammate_distance = TeammateDistanceTo(teammate, ball);
    if (opponent != Unum_Unknown)
        opponent_distance = OpponentDistanceTo(opponent, ball);

    if (teammate_distance < opponent_distance && BallKickableForTeammate(teammate, buffer))
        return teammate;
    else if (opponent_distance < teammate_distance && BallKickableForOpponent(opponent, buffer))
        return -opponent;
    else
        return Unum_Unknown;
}

/*****************************************************************************************/

Unum PositionInfo::TeammateWithBall(float buffer)
{
    Unum player = PlayerWithBall(buffer);
    if (player > 0)
        return player;
    else
        return Unum_Unknown;
}

/*****************************************************************************************/

Unum PositionInfo::OpponentWithBall(float buffer)
{
    Unum player = PlayerWithBall(buffer);
    if (player < 0)
        return -player;
    else
        return Unum_Unknown;
}

/*****************************************************************************************/

char PositionInfo::TeamWithBall(float buffer)
{
    Unum player = PlayerWithBall(buffer);
    if (player == Unum_Unknown)
        return '?';
    if (player > 0)
        return MySide;
    if (player < 0)
        return TheirSide;
    my_error("player with ball needs to say something");
    return '*';
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

int PositionInfo::PlayerPredictedCyclesToPoint(char side, Unum num, Vector pt,
                                               float dash_power, float buffer)
{
    if (side == MySide && num == MyNumber)
        return PredictedCyclesToPoint(pt, dash_power);

    if (!PlayerPositionValid(side, num))
    {
        my_error("Can't predict cycles to point if position invalid %c %d", side, num);
        return 10000;
    }

    LogAction7(210, "PlayerPredCyclesToPoint: %d, (%.1f, %.1f) %.1f, %.1f",
               (side == MySide ? num : -num), pt.x, pt.y, dash_power, buffer);

    float bodyang;
    if (PlayerBodyAngleValid(side, num))
        bodyang = PlayerAbsoluteBodyAngle(side, num);
    else
        bodyang = (pt - PlayerAbsolutePosition(side, num)).dir();
    Vector position = PlayerAbsolutePosition(side, num);
    Vector velocity;
    if (!PlayerVelocityValid(side, num))
        velocity = 0;
    else
        velocity = MyVel();

    for (int i = 0; TRUE; i++)
    {
        if (position.dist(pt) <= buffer)
        {
            LogAction6(210, "PlayerPredicatedPosition: %d at (%.1f, %.1f) spd: %.2f. I'm there",
                       i, position.x, position.y, velocity.mod());
            return i;
        }

        /* decide if they should turn */
        float targ_ang = (pt - position).dir() - bodyang;
        NormalizeAngleDeg(&targ_ang);
        if (fabs(targ_ang) > CP_max_go_to_point_angle_err)
        {
            /* turning */
            float this_turn = MinMax(-EffectiveTurn(SP_max_moment, velocity.mod()),
                                     targ_ang,
                                     EffectiveTurn(SP_max_moment, velocity.mod()));
            bodyang += this_turn;
            LogAction7(210, "PlayerPredicatedPosition: %d at (%.1f, %.1f) spd: %.2f. turning %.2f",
                       i, position.x, position.y, velocity.mod(), this_turn);
        }
        else
        {
            /* dashing */
            velocity += Polar2Vector(dash_power * SP_dash_power_rate, bodyang);
            LogAction6(210, "PlayerPredicatedPosition: %d at (%.1f, %.1f) spd: %.2f. dashing",
                       i, position.x, position.y, velocity.mod());
        }

        if (velocity.mod() > SP_player_speed_max)
            velocity *= (SP_player_speed_max / velocity.mod());

        position += velocity;
        velocity *= SP_player_decay;
    }
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

PlayerObject *PositionInfo::GetTeammate(Unum num)
{
    if (num == Unum_Unknown)
        my_error("Shouldn't get a teammate with unknown num");
    if (num == MyNumber)
        my_error("Shouldn't get self from MyPlayer");

    for (int i = 0; i < num_my_players; i++)
        if (MyPlayer[i]->unum == num)
            return MyPlayer[i];

    return NULL;
}

/*********************************************************************************/

PlayerObject *PositionInfo::GetOpponent(Unum num)
{
    if (num == Unum_Unknown)
        my_error("Shouldn't get an opponent with unknown num");

    for (int i = 0; i < num_their_players; i++)
        if (TheirPlayer[i]->unum == num)
            return TheirPlayer[i];

    return NULL;
}

/*********************************************************************************/

PlayerObject *PositionInfo::GetPlayer(char side, Unum num)
{
    if (side == MySide)
        return GetTeammate(num);
    else if (side == TheirSide)
        return GetOpponent(num);
    else
        my_error("Can't get a player from an unknown side");

    return NULL;
}

/*********************************************************************************/

PlayerObject *PositionInfo::GetNewPlayer(char side, Unum num)
{
    while (num_free_players <= 0 ||
           (side == MySide && num_my_players >= SP_team_size - 1) ||
           (side == TheirSide && num_their_players >= SP_team_size))
    {
        if (!ForgetAPlayer(side)) /* All the players on that side are valid */
            return NULL;
    }

    PlayerObject *player = FreePlayer[--num_free_players];

    if (side == MySide)
        MyPlayer[num_my_players++] = player;
    else if (side == TheirSide)
        TheirPlayer[num_their_players++] = player;
    else /* side == '?' */
        TeamlessPlayer[num_teamless_players++] = player;

    player->side = side;
    player->unum = num;

    return player;
}

/*********************************************************************************/

Bool PositionInfo::ForgetAPlayer(char side)
{
    /* char original_side = side; */

    /* If there aren't enough on the team, forget a teamless player */
    if (side == MySide && num_my_players < SP_team_size - 1)
        side = '?';
    if (side == TheirSide && num_their_players < SP_team_size)
        side = '?';

    /* get rid of the one with unknown numbers that have lowest conf */
    int least_conf_index = -1;
    float least_conf = CP_max_conf + 1; /* Must forget a player, even if all at max conf */
    for (int i = 0; i < num_players; i++)
    {

        if (side == 'f')
            continue;
        if (side != '?' && side != Player[i].side)
            continue;

        if ((Player[i].unum == Unum_Unknown || Player[i].side == '?') && Player[i].pos_valid() < least_conf)
        {
            least_conf = Player[i].pos_valid();
            least_conf_index = i;
        }
    }

    if (least_conf_index >= 0)
    {
        side = Player[least_conf_index].side;
        Player[least_conf_index].reset();
    }
    else if (side == MySide && num_my_players > SP_team_size - 1)
    {
        if (ResetMyDuplicatePlayers() == FALSE)
        {
            my_error("There should be a duplicate MyPlayer");
            return FALSE; /* else do the clean below */
        }
    }
    else if (side == TheirSide && num_their_players > SP_team_size)
    {
        if (ResetTheirDuplicatePlayers() == FALSE)
        {
            my_error("There should be a duplicate TheirPlayer");
            return FALSE; /* else do the clean below */
        }
    }
    else
        return FALSE;

    CleanAllPlayers();
    /* if      ( side == MySide    ) CleanMyPlayers();
    else if ( side == TheirSide ) CleanTheirPlayers();
    else if ( side == '?'       ) CleanTeamlessPlayers();
    else my_error("Which side?"); */

    return TRUE;
}

/*********************************************************************************/

void PositionInfo::CleanMyPlayers()
{
    int new_num_my_players = 0;

    for (int i = 0; i < num_my_players; i++)
    {
        if (MyPlayer[i]->pos_valid())
            MyPlayer[new_num_my_players++] = MyPlayer[i];
        else
        {
            MyPlayer[i]->side = 'f';
            MyPlayer[i]->unum = Unum_Unknown;
            MyPlayer[i]->reset();
            FreePlayer[num_free_players++] = MyPlayer[i];
        }
    }

    num_my_players = new_num_my_players;

    while (num_my_players > SP_team_size - 1)
    {
        /* my_stamp; printf("%d of my players\n",num_my_players);    */
        if (ForgetAPlayer(MySide) == FALSE)
        { /* recurses through CleanMyPlayers */
            my_error("Should be able to forget a teammate");
            printf("Number of players (%d %d %d %d)\n",
                   num_my_players, num_their_players, num_teamless_players, num_free_players);
            for (int i = 0; i < num_players; i++)
                printf("%c %d   %f\n", Player[i].side, Player[i].unum, Player[i].pos_valid());
        }
    }

    if (num_my_players > SP_team_size - 1)
    {
        my_error("Too many of my players %d", num_my_players);
        /* my_stamp;
        for (int i=0; i<num_my_players; i++)
          printf("%d %.1f  ",MyPlayer[i]->unum,MyPlayer[i]->pos_valid() );
        printf("\n");
        dump_core("dump");*/
    }
}

/*********************************************************************************/

void PositionInfo::CleanTheirPlayers()
{
    int new_num_their_players = 0;

    for (int i = 0; i < num_their_players; i++)
    {
        if (TheirPlayer[i]->pos_valid())
            TheirPlayer[new_num_their_players++] = TheirPlayer[i];
        else
        {
            TheirPlayer[i]->side = 'f';
            TheirPlayer[i]->unum = Unum_Unknown;
            TheirPlayer[i]->reset();
            FreePlayer[num_free_players++] = TheirPlayer[i];
        }
    }
    num_their_players = new_num_their_players;

    while (num_their_players > SP_team_size)
    {
        /* my_stamp; printf("%d of their players\n",num_their_players);    */
        if (ForgetAPlayer(TheirSide) == FALSE)
        { /* recurses through CleanTheirPlayers */
            my_error("Should be able to forget an opponent");
            printf("Number of players (%d %d %d %d)\n",
                   num_my_players, num_their_players, num_teamless_players, num_free_players);
            for (int i = 0; i < num_players; i++)
                printf("%c %d   %f\n", Player[i].side, Player[i].unum, Player[i].pos_valid());
        }
    }

    if (num_their_players > SP_team_size)
    {
        my_error("Too many of their players %d", num_their_players);
        /*my_stamp;
        for (int i=0; i<num_their_players; i++)
          printf("%d %.1f  ",TheirPlayer[i]->unum,TheirPlayer[i]->pos_valid() );
        printf("\n");
        dump_core("dump");*/
    }
}

/*********************************************************************************/

void PositionInfo::CleanTeamlessPlayers()
{
    int new_num_teamless_players = 0;

    for (int i = 0; i < num_teamless_players; i++)
    {
        if (TeamlessPlayer[i]->pos_valid())
        {
            /* player may have been identified as being on one team or another */
            if (TeamlessPlayer[i]->side == MySide)
            {
                MyPlayer[num_my_players++] = TeamlessPlayer[i];
            }
            else if (TeamlessPlayer[i]->side == TheirSide)
            {
                TheirPlayer[num_their_players++] = TeamlessPlayer[i];
            }
            else if (TeamlessPlayer[i]->side == '?')
                TeamlessPlayer[new_num_teamless_players++] = TeamlessPlayer[i];
            else
                my_error("Teamless players should have side '?'");
        }
        else
        {
            TeamlessPlayer[i]->side = 'f';
            TeamlessPlayer[i]->unum = Unum_Unknown;
            TeamlessPlayer[i]->reset();
            FreePlayer[num_free_players++] = TeamlessPlayer[i];
        }
    }

    num_teamless_players = new_num_teamless_players;

    while (num_teamless_players > num_players)
        if (ForgetAPlayer('?') == FALSE) /* recurses through CleanTeamlessPlayers */
            my_error("Should be able to forget a teamless player");

    if (num_teamless_players > num_players)
    {
        my_error("Too many of teamless players %d", num_teamless_players);
        my_stamp;
        for (int i = 0; i < num_teamless_players; i++)
            printf("%d %.1f  ", TeamlessPlayer[i]->unum, TeamlessPlayer[i]->pos_valid());
        printf("\n");
    }
}

/*********************************************************************************/

void PositionInfo::CleanAllPlayers()
{
    CleanTeamlessPlayers(); /* has to be first to move over reconciled players */
    CleanMyPlayers();
    CleanTheirPlayers();
}

/*********************************************************************************/

Bool PositionInfo::ResetMyDuplicatePlayers()
{
    if (num_my_players <= SP_team_size - 1)
        my_error("Shouldn't be duplicates");

    int MyPlayerIndex[MAX_PLAYERS]; /* Map uniform number to index in MyPlayer */

    for (int num = 1; num <= SP_team_size; num++)
        MyPlayerIndex[num] = -1;

    Bool DuplicateFound = FALSE;
    Unum unum;
    for (int i = 0; i < num_my_players; i++)
    {
        unum = MyPlayer[i]->unum;
        if (unum == Unum_Unknown)
            my_error("Catch unknowns before duplicates (my)");

        if (unum == MyNumber)
        {
            MyPlayer[i]->reset();
            DuplicateFound = TRUE;
            my_error("How did my number get in the mix?", unum);
        }
        else if (MyPlayerIndex[unum] == -1)
            MyPlayerIndex[unum] = i;
        else if (MyPlayer[MyPlayerIndex[unum]]->pos_valid() >= MyPlayer[i]->pos_valid())
        {
            MyPlayer[i]->reset();
            DuplicateFound = TRUE;
            my_error("There were 2 %d's on my side", unum);
        }
        else
        {
            MyPlayer[MyPlayerIndex[unum]]->reset();
            MyPlayerIndex[unum] = i;
            DuplicateFound = TRUE;
            my_error("There were 2 %d's on my side", unum);
        }
    }

    return DuplicateFound;
}

/*********************************************************************************/

Bool PositionInfo::ResetTheirDuplicatePlayers()
{
    if (num_their_players <= SP_team_size)
        my_error("Shouldn't be duplicates");

    int TheirPlayerIndex[MAX_PLAYERS]; /* Map uniform number to index in TheirPlayer */

    for (int num = 1; num <= SP_team_size; num++)
        TheirPlayerIndex[num] = -1;

    Bool DuplicateFound = FALSE;
    Unum unum;
    for (int i = 0; i < num_their_players; i++)
    {
        unum = TheirPlayer[i]->unum;
        if (unum == Unum_Unknown)
            my_error("Catch unknowns before duplicates (their)");

        if (TheirPlayerIndex[unum] == -1)
            TheirPlayerIndex[unum] = i;
        else if (TheirPlayer[TheirPlayerIndex[unum]]->pos_valid() >= TheirPlayer[i]->pos_valid())
        {
            TheirPlayer[i]->reset();
            DuplicateFound = TRUE;
            my_error("There were 2 %d's on their side", unum);
        }
        else
        {
            TheirPlayer[TheirPlayerIndex[unum]]->reset();
            TheirPlayerIndex[unum] = i;
            DuplicateFound = TRUE;
            my_error("There were 2 %d's on their side", unum);
        }
    }

    return DuplicateFound;
}

/*********************************************************************************/

void PositionInfo::ClearSeenInfo()
{
    /* In case we need to forget a sight because a new one came in before we could update */
    Ball.clear_seen();
    for (int i = 0; i < num_players; i++)
        Player[i].clear_seen();
    num_unknown_players = 0;

    for (int i = 0; i < num_seen_markers; i++)
        Marker[SeenMarker[i]].clear_seen(); /* Not necessarily needed... */
    num_seen_markers = 0;
    ClosestMarker = ClosestMotionMarker = No_Marker;
    SeenLine = SL_No_Line;
}

/*********************************************************************************/

MarkerType PositionInfo::ClosestGoal()
{
    if (MyConf() && MyX() < 0)
        return RM_My_Goal;
    else
        return RM_Their_Goal;
}

/*********************************************************************************/

MarkerType PositionInfo::ClosestFlagTo()
{
    return No_Marker;
}

/*********************************************************************************/

PlayerObject *PositionInfo::ClosestPlayerObjectTo(Vector gpos)
{
    /* Need to check teammates, opponents, teamless */
    /* return NULL if no-one's within player-max-speed (sqr) of the right place */

    int i;
    PlayerObject *player = NULL;
    float min_dist_sqr = 40000; /* More than pitch_length*pitch_length + pitch_width*pitch_width */
    /* If player is farther than 4 times the maximum single-cycle movement, it can't be the same one */
    float max_dist_sqr = (SP_player_speed_max * CP_max_player_move_factor) *
                         (SP_player_speed_max * CP_max_player_move_factor);
    float dist_sqr;

    for (i = 0; i < num_my_players; i++)
    {
        /* not in memory yet, so not valid -- must be different player */
        if (!(MyPlayer[i]->pos_valid()))
            continue;

        dist_sqr = gpos.dist2(MyPlayer[i]->get_abs_pos());
        if (dist_sqr < max_dist_sqr && dist_sqr < min_dist_sqr)
        {
            min_dist_sqr = dist_sqr;
            player = MyPlayer[i];
        }
    }
    for (i = 0; i < num_their_players; i++)
    {
        /* not in memory yet, so not valid -- must be different player */
        if (!(TheirPlayer[i]->pos_valid()))
            continue;

        dist_sqr = gpos.dist2(TheirPlayer[i]->get_abs_pos());
        if (dist_sqr < max_dist_sqr && dist_sqr < min_dist_sqr)
        {
            min_dist_sqr = dist_sqr;
            player = TheirPlayer[i];
        }
    }
    for (i = 0; i < num_teamless_players; i++)
    {
        /* not in memory yet, so not valid -- must be different player */
        if (!(TeamlessPlayer[i]->pos_valid()))
            continue;

        dist_sqr = gpos.dist2(TeamlessPlayer[i]->get_abs_pos());
        if (dist_sqr < max_dist_sqr && dist_sqr < min_dist_sqr)
        {
            min_dist_sqr = dist_sqr;
            player = TeamlessPlayer[i];
        }
    }

    return player;
}

/*********************************************************************************/

PlayerObject *PositionInfo::ClosestTeammateObjectTo(Vector gpos)
{
    /* Need to check teammates, teamless */
    /* return NULL if no-one's within player-max-speed (sqr) of the right place */

    int i;
    PlayerObject *player = NULL;
    float min_dist_sqr = 40000; /* More than pitch_length*pitch_length + pitch_width*pitch_width */
    /* If player is farther than 4 times the maximum single-cycle movement, it can't be the same one */
    float max_dist_sqr = (SP_player_speed_max * CP_max_player_move_factor) *
                         (SP_player_speed_max * CP_max_player_move_factor);
    float dist_sqr;

    for (i = 0; i < num_my_players; i++)
    {
        /* not in memory yet, so not valid -- must be different player */
        if (!(MyPlayer[i]->pos_valid()))
            continue;

        dist_sqr = gpos.dist2(MyPlayer[i]->get_abs_pos());
        if (dist_sqr < max_dist_sqr && dist_sqr < min_dist_sqr)
        {
            min_dist_sqr = dist_sqr;
            player = MyPlayer[i];
        }
    }
    for (i = 0; i < num_teamless_players; i++)
    {
        /* not in memory yet, so not valid -- must be different player */
        if (!(TeamlessPlayer[i]->pos_valid()))
            continue;

        dist_sqr = gpos.dist2(TeamlessPlayer[i]->get_abs_pos());
        if (dist_sqr < max_dist_sqr && dist_sqr < min_dist_sqr)
        {
            min_dist_sqr = dist_sqr;
            player = TeamlessPlayer[i];
        }
    }

    return player;
}

/*********************************************************************************/

PlayerObject *PositionInfo::ClosestOpponentObjectTo(Vector gpos)
{
    /* Need to check opponents, teamless */
    /* return NULL if no-one's within player-max-speed (sqr) of the right place */

    int i;
    PlayerObject *player = NULL;
    float min_dist_sqr = 40000; /* More than pitch_length*pitch_length + pitch_width*pitch_width */
    /* If player is farther than 4 times the maximum single-cycle movement, it can't be the same one */
    float max_dist_sqr = (SP_player_speed_max * CP_max_player_move_factor) *
                         (SP_player_speed_max * CP_max_player_move_factor);
    float dist_sqr;

    for (i = 0; i < num_their_players; i++)
    {
        /* not in memory yet, so not valid -- must be different player */
        if (!(TheirPlayer[i]->pos_valid()))
            continue;

        dist_sqr = gpos.dist2(TheirPlayer[i]->get_abs_pos());
        if (dist_sqr < max_dist_sqr && dist_sqr < min_dist_sqr)
        {
            min_dist_sqr = dist_sqr;
            player = TheirPlayer[i];
        }
    }
    for (i = 0; i < num_teamless_players; i++)
    {
        /* not in memory yet, so not valid -- must be different player */
        if (!(TeamlessPlayer[i]->pos_valid()))
            continue;

        dist_sqr = gpos.dist2(TeamlessPlayer[i]->get_abs_pos());
        if (dist_sqr < max_dist_sqr && dist_sqr < min_dist_sqr)
        {
            min_dist_sqr = dist_sqr;
            player = TeamlessPlayer[i];
        }
    }

    return player;
}

/*********************************************************************************/

PlayerObject *PositionInfo::ClosestPlayerObjectTo(char side, Vector gpos)
{
    if (side == MySide)
        return ClosestTeammateObjectTo(gpos);
    else if (side == TheirSide)
        return ClosestOpponentObjectTo(gpos);
    else /* side == '?' */
        return ClosestPlayerObjectTo(gpos);
}

/*********************************************************************************/

/* SMURF: as a hack, if n is Unum_Teamless, it does it for closest
   teamless player */
Bool PositionInfo::BallKickableForPlayer(char s, Unum n, float buffer)
{
    Vector pos;
    if (n == Unum_Teamless)
    {
        if (NumTeamlessPlayers() < 1)
            my_error("Can't tell kickable for teamless if there aren;t any");
        pos = ClosestTeamlessPlayerPosition();
    }
    else
    {
        if (!PlayerPositionValid(s, n))
            my_error("Can't tell you if they can kick it if we don;t know where they are");
        pos = PlayerAbsolutePosition(s, n);
    }
    return ((pos - BallAbsolutePosition()).mod() <
            SP_kickable_area - buffer)
               ? TRUE
               : FALSE;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

void PositionInfo::update()
{
    if (NewSight && LastSightTime < CurrentTime - 1)
    {
        /* Special case--ignore sights from just before change to play_on_mode
           if they seem to be 2 cycles ago.  Better would be to adjust the sight
           time, but then sight times of all other objects also need to be adjusted */
        NewSight = FALSE;
        if (LastSightTime < CurrentTime - 2)
            my_error("last sight shouldn't be so out of date");
    }

    /* before updating from any sights: */
    my_last_body_ang = MyBodyAng();
    my_last_neck_global_ang = MyNeckGlobalAng();
    my_last_neck_rel_ang = MyNeckRelAng();

    /* see the lengthy comment above titled POSITION BASED VELOCITY ESTIMATION */
    /* this next part used to be just for position based velocity estimation, but now
       we use it for ball velocity invalidation as well
       The comment referred to previously also applies to ball velocity invalidation */
    Vector estimated_pos;
    Bool estimated_pos_valid = FALSE;
    // if (CP_use_new_position_based_vel && NewSight && MyConf()) {
    if (NewSight && MyConf())
    {
        update_self_estimate(LastSightTime);
        estimated_pos = MyPos();
        estimated_pos_valid = TRUE;
        // cout << "Time: " << CurrentTime.t << "\tEstPos: " << estimated_pos << endl;
    }

    /* So my velocity and position match for computing ball, player positions */
    /* Updates my position, angle and velocity at LastSightTime               */
    if (NewSight && LastSightTime < CurrentTime)
        update_self_seen(LastSightTime);
    else
        update_self(CurrentTime);

    /* see the lengthy comment above titled POSITION BASED VELOCITY ESTIMATION */
    /* this next part used to be just for position based velocity estimation, but now
       we use it for ball velocity invalidation as well
       The comment referred to previously also applies to ball velocity invalidation */
    // if (CP_use_new_position_based_vel && NewSight && MyConf() && estimated_pos_valid) {
    if (NewSight && MyConf() && estimated_pos_valid)
    {
        sight_position_correction = MyPos() - estimated_pos;
        sight_position_correction_time = LastSightTime;
    }

    update_ball(CurrentTime);
    update_players(CurrentTime);

    /* Brings me from time prev time to CurrentTime if necessary */
    if (MyUpdateTime() < CurrentTime)
        update_self_estimate(CurrentTime);

    update_stamina(CurrentTime);
    update_offside_lines();

#ifndef NO_ACTION_LOG
    if (CP_save_action_log_level >= 175)
    {
        char outstring[200], *pc;

        if (!Mem->MyConf())
            sprintf(outstring, "My Pos: (??, ??)\tangle: ??\trel_neck: %.2f",
                    MyNeckRelAng());
        else
            sprintf(outstring, "My Pos: (%.2f, %.2f)\tangle: %.2f\trel_neck: %.2f\tvel: (%.2f %.2f)\tconf: %.2f\tstamina: %.2f",
                    MyX(), MyY(), MyBodyAng(), MyNeckRelAng(), MySpeed(), MyDir(), MyConf(), MyStamina());
        LogAction2(175, outstring);

        if (!Mem->BallPositionValid())
            sprintf(outstring, "Ball Pos: (??, ??)");
        else if (!Mem->BallVelocityValid())
            sprintf(outstring, "Ball Pos: (%.2f, %.2f)conf: %.2f\tVel: (??, ??)",
                    BallX(), BallY(), BallPositionValid());
        else
            sprintf(outstring, "Ball Pos: (%.2f, %.2f)conf: %.2f\tVel: (%.2f,%.2f)dx/dy = (%.2f,%.2f)sp/head conf: %.2f",
                    BallX(), BallY(), BallPositionValid(),
                    BallAbsoluteVelocity().x, BallAbsoluteVelocity().y,
                    BallSpeed(), BallAbsoluteHeading(),
                    BallVelocityValid());
        LogAction2(175, outstring);

        strcpy(outstring, "Player pos known: team: ");
        pc = outstring + strlen(outstring);
        for (int num = 1; num <= SP_team_size; num++)
            *(pc++) = (TeammatePositionValid(num) ? char_for_num(num) : '_');
        *pc = 0; /* null terminate */
        strcat(outstring, "\topp: ");
        pc = outstring + strlen(outstring);
        for (int num = 1; num <= SP_team_size; num++)
            *(pc++) = (OpponentPositionValid(num) ? char_for_num(num) : '_');
        *pc = 0; /* null terminate */
        LogAction2(200, outstring);

        if (NewSight)
        {
            /* record what we saw for this new sight */
            sprintf(outstring, "Sight at %d.%d: %c%c team: ", LastSightTime.t,
                    LastSightTime.s,
                    (Ball.GetSeenTime() == LastSightTime ? 'B' : '_'),
                    (Ball.GetSeenMovingTime() == LastSightTime ? 'v' : '_'));
            pc = outstring + strlen(outstring);
            for (int num = 1; num <= SP_team_size; num++)
            {
                if (num != MyNumber && GetTeammate(num))
                    *(pc++) = (GetTeammate(num)->GetSeenTime() == LastSightTime ? char_for_num(num) : '_');
                else
                    *(pc++) = '_';
            }
            *pc = 0; /* null terminate */
            strcat(outstring, "\topp: ");
            pc = outstring + strlen(outstring);
            for (int num = 1; num <= SP_team_size; num++)
            {
                if (GetOpponent(num))
                    *(pc++) = (GetOpponent(num)->GetSeenTime() == LastSightTime ? char_for_num(num) : '_');
                else
                    *(pc++) = '_';
            }
            *pc = 0; /* null terminate */
            LogAction2(175, outstring);
        }
    }
#endif

    NewAction = FALSE;
    NewSight = FALSE;
}

/*********************************************************************************/

void PositionInfo::update_self_seen(Time time)
{
    /* Brings me to the best known position and speed as of time x through vision and
       past values */
    if (!NewSight || time < CurrentTime - 1)
        my_error("No new sight with which to update %d %d (%d %d)",
                 time.t, time.s, LastStartClockTime.t, LastStartClockTime.s);

    update_self_neck_rel_ang(time);

    if (SeenLine != SL_No_Line && ClosestMarker != No_Marker)
    {
        Fieldline[SeenLine].sanitize_times();
        Marker[ClosestMarker].sanitize_times();

        SetMyBodyAng(GetNormalizeAngleDeg(Fieldline[SeenLine].get_my_neck_global_ang() - MyNeckRelAng()));
        SetMyPos(Marker[ClosestMarker].get_my_pos(MyNeckGlobalAng()), time);

        if ((MyVelConf() < CP_max_conf || my_vel_time != time) && SensedInfoKnown(time))
            SetMyVel(Polar2Vector(GetMySensedSpeed(time), MyBodyAng()), time);
    }
    // else my_error("NO MARKER");  /* not really an error---can estimate */

    for (int i = 0; i < num_seen_markers; i++)
        Marker[SeenMarker[i]].reset(); /* Not necessarily needed... */
    num_seen_markers = 0;
    ClosestMarker = ClosestMotionMarker = No_Marker;
    SeenLine = SL_No_Line;
}

/*********************************************************************************/

void PositionInfo::update_self(Time time)
{
    if (NewSight)
    {
        if (LastSightTime != CurrentTime)
            my_error("shouldn't be here");
        update_self_seen(time);
    }

    /* Brings unknown values up to date */
    if (MyUpdateTime() < time)
        update_self_estimate(time);
}

/*********************************************************************************/

void PositionInfo::update_ball(Time time)
{
    Ball.update(time);
}

/*********************************************************************************/

void PositionInfo::update_players(Time time)
{
    int i;

    /* Assume my position's already updated here */
    if (NewSight && MyConf())
        reconcile_unknown_players();

    for (i = 0; i < num_my_players; i++)
        MyPlayer[i]->update(time);
    for (i = 0; i < num_their_players; i++)
        TheirPlayer[i]->update(time);
    for (i = 0; i < num_teamless_players; i++)
        TeamlessPlayer[i]->update(time);

    CleanAllPlayers();
    num_unknown_players = 0;

    if (num_free_players + num_my_players + num_their_players + num_teamless_players != num_players)
        my_error("Number of players doesn't add up (%d %d %d %d)",
                 num_my_players, num_their_players, num_teamless_players, num_free_players);
}

/*********************************************************************************/

void PositionInfo::reconcile_unknown_players()
{
    if (!MyConf() || MyPosTime() != LastSightTime)
        my_error("Can't reconcile unknown players if not localized");

    PlayerObject *player;
    TempPlayerObject *unknown_player;

    char s;
    float d;
    AngleDeg a;
    Time t;

    for (int i = 0; i < num_unknown_players; i++)
    {

        unknown_player = &(UnknownPlayer[i]);
        s = unknown_player->side;
        d = unknown_player->dist;
        a = unknown_player->ang_from_neck;
        t = unknown_player->time;

        if (t != LastSightTime)
            continue; /* For some reason, some old players get to here */

        Vector rel_pos = Polar2Vector(d, a);
        Vector global_pos = MyPos() + rel_pos.rotate(MyNeckGlobalAng());

        if ((player = ClosestPlayerObjectTo(s, global_pos)) == NULL)
            player = GetNewPlayer(s, Unum_Unknown);
        else if (s != '?' && player->side == '?')
            player->side = s; /* know the teamless player's side now  */
                              /* need to call CleanTeamlessPlayers() before doing much else
                                 (like CleanMyPlayers(), ForgetAPlayer(), etc.*/
        if (player != NULL && (player->unum != MyNumber || player->side != MySide))
            player->set_polar_from_neck(d, a, t);
    }
}

/*********************************************************************************/

void PositionInfo::update_offside_lines()
{
    if (!SP_use_offside)
        return;

    float first = 0, second = 0, tmp;

    for (int i = 0; i < num_my_players; i++)
    {
        if (MyPlayer[i]->pos_valid() && MyPlayer[i]->get_abs_pos().Behind(second))
        {
            second = MyPlayer[i]->get_abs_pos().x;
            if (second < first)
            {
                tmp = first;
                first = second;
                second = tmp;
            }
        }
    }
    if (BallPositionValid() && BallAbsolutePosition().Behind(second))
        their_offside_line = BallAbsolutePosition().x;
    else
        their_offside_line = second;

    first = second = 0;

    for (int i = 0; i < num_their_players; i++)
    {
        if (TheirPlayer[i]->pos_valid() && TheirPlayer[i]->get_abs_pos().InFrontOf(second))
        {
            second = TheirPlayer[i]->get_abs_pos().x;
            if (second > first)
            {
                tmp = first;
                first = second;
                second = tmp;
            }
        }
    }
    if (BallPositionValid() && BallAbsolutePosition().InFrontOf(second))
        my_offside_line = BallAbsolutePosition().x;
    else
        my_offside_line = second;
}

/*********************************************************************************/

Bool PositionInfo::OffsidePosition(float x, char side)
{
    if (!SP_use_offside)
        return FALSE;

    if (side == MySide)
        return (x > my_offside_line) ? TRUE : FALSE;
    else if (side == TheirSide)
        return (x < their_offside_line) ? TRUE : FALSE;
    else
        my_error("Can't tell offside if don't know the player's team");
    return FALSE;
}

/*********************************************************************************/

Bool PositionInfo::TeammateInOffsidePosition(Unum num)
{
    if (!SP_use_offside)
        return FALSE;
    if (!TeammatePositionValid(num))
        my_error("Can't tell if teammate offside--not valid");

    return OffsidePosition(TeammateAbsolutePosition(num), MySide);
}

/*********************************************************************************/

Bool PositionInfo::OpponentInOffsidePosition(Unum num)
{
    if (!SP_use_offside)
        return FALSE;
    if (!OpponentPositionValid(num))
        my_error("Can't tell if opponent offside--not valid");

    return OffsidePosition(OpponentAbsolutePosition(num), TheirSide);
}

/*********************************************************************************/

Bool PositionInfo::PlayerInOffsidePosition(char side, Unum num)
{
    if (!SP_use_offside)
        return FALSE;

    if (side == MySide)
        return TeammateInOffsidePosition(num);
    else if (side == TheirSide)
        return OpponentInOffsidePosition(num);
    else
        my_error("Can't tell offside if don't know the player's team");
    return FALSE;
}

/*********************************************************************************/

Unum PositionInfo::TeammateOffsideIfIKick()
{
    if (!SP_use_offside)
        return Unum_Unknown;

    for (int teammate = 1; teammate <= SP_team_size; teammate++)
    {
        if (teammate != MyNumber && TeammatePositionValid(teammate) &&
            TeammateDistance(teammate) < SP_offside_area && TeammateInOffsidePosition(teammate))
            return teammate;
    }

    return Unum_Unknown;
}

/*****************************************************************************************/

float PositionInfo::XToAdjustForOffsideX(float x, float buffer)
{
    if (!SP_use_offside)
        return x;

    float back_x = -SP_pitch_length / 2;
    if (PullOffside)
        back_x = PullOffsidePosition + buffer;

    float front_x = my_offside_line;

    if (BallPositionValid())
        back_x = Min(back_x, BallX());

    if (x < back_x) /* bring players up */
        x = back_x;
    else if (x <= my_offside_line)
        /* move midfielders back in proportion to difference between pull_x and my_offside_line
           relative to the length of the field */
        x = ((x + SP_pitch_length / 2) / SP_pitch_length) * (front_x - back_x) + back_x;

    if (OffsidePosition(x, MySide))
        x = XToOnsideX(x, buffer);

    return x;
}

/*****************************************************************************************/

Rectangle PositionInfo::RectangleToAdjustForOffsideRectangle(Rectangle *rect, float buffer)
{
    if (!SP_use_offside)
        return *rect;

    float right = XToAdjustForOffsideX(rect->RightX(), buffer);
    float left = Min(right - 5, XToAdjustForOffsideX(rect->LeftX(), buffer));

    float top = rect->TopY();
    float bottom = rect->BottomY();

    return Rectangle(left, right, top, bottom);
}

/*****************************************************************************************/

Vector PositionInfo::PositionToPullOffsidePosition(Vector pos, float buffer)
{
    if (!SP_use_offside)
        return pos;

    float pull_x = PullOffsidePosition + buffer;
    if (BallPositionValid())
        pull_x = Min(pull_x, BallX());

    if (pos.x < pull_x) /* bring players up */
        pos.x = pull_x;
    else if (!OffsidePosition(pos, MySide))
        /* move midfielders back in proportion to difference between pull_x and my_offside_line
           relative to the length of the field */
        pos.x = ((pos.x + SP_pitch_length / 2) / (my_offside_line + SP_pitch_length / 2)) *
                    (my_offside_line - pull_x) +
                pull_x;

    return pos;
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

void PositionInfo::VerifyDash(float *dash_power)
{
    switch (PlayMode)
    {
    case PM_Their_Kick_Off:
    case PM_Their_Kick_In:
    case PM_Their_Free_Kick:
    case PM_Their_Offside_Kick:
    case PM_Their_Corner_Kick: /* Don't waste stamina trying to get closer to ball */
        if (MyConf() && BallPositionValid() &&
            (MyPos() + NewVelFromDash(MyVel(), *dash_power)).dist(BallAbsolutePosition()) < SP_free_kick_buffer)
            *dash_power = 0;
        break;
    case PM_Their_Goal_Kick: /* Don't waste stamina trying to get closer to ball */
        if (MyConf() && BallPositionValid() &&
            TheirPenaltyArea.IsWithin(MyPos() + NewVelFromDash(MyVel(), *dash_power)))
            *dash_power = 0;
        break;
    default:;
    }

    PlayerInfo::VerifyDash(dash_power);
}

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

int PositionInfo::SortPlayersBy(char side, char KeyFunc, float KeyNum, Unum *players)
{
    int result = 0; /*Number of players sorted */

    float (PositionInfo::*KeyFunction)(char, Unum);
    /* taking player angles from body */
    KeyFunction =
        ((KeyFunc == 'd') ? &PositionInfo::PlayerDistance : &PositionInfo::PlayerAngleFromBody);

    int num = ((side == 'b') ? SP_team_size * 2 : SP_team_size); /* Make aux array big
                                        enough */
    float *vals;
    vals = new float[num];

    char team = ((side == 't') ? TheirSide : MySide);
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (PlayerPositionValid(team, i))
        {
            players[result] = i;
            vals[result] = fabs(KeyNum - (this->*KeyFunction)(team, i)); /* use diff from key*/
            result++;
        }
    }

    if (side == 'b')
    { /* Need to put in Their  team too */
        team = TheirSide;
        for (int i = 1; i <= SP_team_size; i++)
        {
            if (PlayerPositionValid(team, i))
            {
                players[result] = i + SP_team_size;                          /* to distinguish from my team */
                vals[result] = fabs(KeyNum - (this->*KeyFunction)(team, i)); /* use diff from key*/
                result++;
            }
        }
    }

    /* Now should have all values in question in vals, with uniform number in
       corresponding position of players ( +TEAM_SIZE for their team if
       side == 'b'):  Just sort em */

    BubbleSort(result, players, vals);
    delete[] vals;
    return result;
}

/*********************************************************************************/

int PositionInfo::SortPlayersByDistanceToPoint(char side, Vector point, Unum *players)
{
    int result = 0; /*Number of players sorted */

    int num = ((side == 'b') ? SP_team_size * 2 : SP_team_size); /* Make aux array big
                                        enough */
    float *vals;
    vals = new float[num];

    char team = ((side == 't') ? TheirSide : MySide);
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (PlayerPositionValid(team, i))
        {
            players[result] = i;
            vals[result] = PlayerDistance2To(team, i, point);
            result++;
        }
    }

    if (side == 'b')
    { /* Need to put in Their  team too */
        team = TheirSide;
        for (int i = 1; i <= SP_team_size; i++)
        {
            if (PlayerPositionValid(team, i))
            {
                players[result] = i + SP_team_size; /* to distinguish from my team */
                vals[result] = PlayerDistance2To(team, i, point);
                result++;
            }
        }
    }

    /* Now should have all values in question in vals, with uniform number in
       corresponding position of players ( +TEAM_SIZE for their team if
       side == 'b'):  Just sort em */

    BubbleSort(result, players, vals);
    delete[] vals;
    return result;
}

/*********************************************************************************/

int PositionInfo::SortPlayersByDistanceToLine(char side, Line line, Unum *players, Bool TestEndPoints, Vector ep1, Vector ep2)
{
    int result = 0; /*Number of players sorted */

    int num = ((side == 'b') ? SP_team_size * 2 : SP_team_size); /* Make aux array big
                                        enough */
    float *vals;
    vals = new float[num];

    Line perp1, perp2;
    float ep_dist2, d1, d2;
    if (TestEndPoints == TRUE)
    {
        perp1 = line.perpendicular(ep1);
        perp2 = line.perpendicular(ep2);
        ep_dist2 = ep1.dist2(ep2);
    }

    char team = ((side == 't') ? TheirSide : MySide);
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (PlayerPositionValid(team, i))
        {
            d1 = PlayerDistance2ToLine(team, i, perp1);
            d2 = PlayerDistance2ToLine(team, i, perp2);
            if (TestEndPoints == TRUE)
            {
                if (d1 > ep_dist2 || d2 > ep_dist2) /* not between points */
                    continue;
            }
            players[result] = i;
            vals[result] = PlayerDistance2ToLine(team, i, line);
            result++;
        }
    }

    if (side == 'b')
    { /* Need to put in Their  team too */
        team = TheirSide;
        for (int i = 1; i <= SP_team_size; i++)
        {
            if (PlayerPositionValid(team, i))
            {
                if (TestEndPoints)
                {
                    if (PlayerDistance2ToLine(team, i, perp1) > ep_dist2 ||
                        PlayerDistance2ToLine(team, i, perp2) > ep_dist2) /* not between points */
                        continue;
                }
                players[result] = i + SP_team_size; /* to distinguish from my team */
                vals[result] = PlayerDistance2ToLine(team, i, line);
                result++;
            }
        }
    }

    /* Now should have all values in question in vals, with uniform number in
       corresponding position of players ( +TEAM_SIZE for their team if
       side == 'b'):  Just sort em */

    BubbleSort(result, players, vals);
    delete[] vals;
    return result;
}

/*********************************************************************************/

int PositionInfo::NumTeammatesWithin(float Dist, Vector of_pos)
{
    int result = 0;
    float Dist2 = Dist * Dist;
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (TeammatePositionValid(i) && TeammateAbsolutePosition(i).dist2(of_pos) <= Dist2 + FLOAT_EPS)
            result++;
    }
    return result;
}

/*********************************************************************************/

int PositionInfo::NumOpponentsWithin(float Dist, Vector of_pos)
{
    int result = 0;
    float Dist2 = Dist * Dist;
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (OpponentPositionValid(i) && OpponentAbsolutePosition(i).dist2(of_pos) <= Dist2 + FLOAT_EPS)
            result++;
    }
    return result;
}

/*********************************************************************************/

int PositionInfo::NumTeammatesWithin(float Dist, AngleDeg Ang, float ofDist, AngleDeg ofAng)
{
    /* ofAng relative to body */
    int result = 0;
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (i != MyNumber && TeammatePositionValid(i) &&
            fabs(TeammateDistance(i) - ofDist) <= Dist + FLOAT_EPS &&
            fabs(GetNormalizeAngleDeg(TeammateAngleFromBody(i) - ofAng)) <= Ang + FLOAT_EPS)
            result++;
    }
    return result;
}

/*********************************************************************************/

int PositionInfo::NumOpponentsWithin(float Dist, AngleDeg Ang, float ofDist, AngleDeg ofAng)
{
    /* ofAng relative to body */
    int result = 0;
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (OpponentPositionValid(i) &&
            fabs(OpponentDistance(i) - ofDist) <= Dist + FLOAT_EPS &&
            fabs(GetNormalizeAngleDeg(OpponentAngleFromBody(i) - ofAng)) <= Ang + FLOAT_EPS)
            result++;
    }
    return result;
}

/*********************************************************************************/

PlayerObject *PositionInfo::GetPlayerWithin(float Dist, Vector ofPos)
{
    float Dist2 = Dist * Dist;
    for (int i = 0; i < num_players; i++)
    {
        if (!(Player[i].unum == MyNumber && Player[i].side == MySide) &&
            Player[i].pos_valid() &&
            Player[i].get_abs_pos().dist2(ofPos) <= Dist2 + FLOAT_EPS)
            return &Player[i];
    }
    return NULL;
}

/*********************************************************************************/

PlayerObject *PositionInfo::GetPlayerWithin(float Dist, AngleDeg Ang, float ofDist, AngleDeg ofAng)
{
    /* ofAng relative to body */
    for (int i = 0; i < num_players; i++)
    {
        if (!(Player[i].unum == MyNumber && Player[i].side == MySide) &&
            Player[i].pos_valid() &&
            fabs(Player[i].get_dist() - ofDist) <= Dist + FLOAT_EPS &&
            fabs(GetNormalizeAngleDeg(Player[i].get_ang_from_body() - ofAng)) <= Ang + FLOAT_EPS)
            return &Player[i];
    }
    return NULL;
}

/*********************************************************************************/

/* end is the end point of the center line of the cone
   wid_dist ratio defines the wid of the cone = wid at dist 1 */
int PositionInfo::NumOpponentsInCone(float wid_dist_ratio, Vector end, Vector vert)
{
    int count = 0;
    Line l = LineFromTwoPoints(vert, end);
    for (Unum opp = 1; opp <= SP_team_size; opp++)
    {
        if (!OpponentPositionValid(opp))
            continue;
        Vector pt = l.ProjectPoint(OpponentAbsolutePosition(opp));
        if (pt.dist2(OpponentAbsolutePosition(opp)) < pt.dist2(vert) * wid_dist_ratio * wid_dist_ratio && l.InBetween(pt, vert, end))
        {
            count++;
        }
    }
    return count;
}

/*********************************************************************************/

/* end is the end point of the center line of the cone
   wid_dist ratio defines the wid of the cone = wid at dist 1 */
int PositionInfo::NumTeammatesInCone(float wid_dist_ratio, Vector end,
                                     Vector vert, Bool IncludeMe)
{
    int count = 0;
    Line l = LineFromTwoPoints(vert, end);
    for (Unum num = 1; num <= SP_team_size; num++)
    {
        if (IncludeMe && num == MyNumber)
            continue;
        if (!TeammatePositionValid(num))
            continue;
        Vector pt = l.ProjectPoint(TeammateAbsolutePosition(num));
        if (pt.dist(TeammateAbsolutePosition(num)) < pt.dist(vert) * wid_dist_ratio * wid_dist_ratio && l.InBetween(pt, vert, end))
        {
            count++;
        }
    }
    return count;
}

/*********************************************************************************/

int PositionInfo::NumPlayersInConeToPlayer(char which,
                                           float wid_dist_ratio, char side,
                                           Unum num, float extra_len, Vector vert)
{
    if (!PlayerPositionValid(side, num))
        my_error("Can't do cone calc to player if we don;t know where he is");

    Vector center = PlayerAbsolutePosition(side, num) - vert;
    float cent_mod = center.mod();
    center *= (cent_mod + extra_len) / cent_mod;
    switch (which)
    {
    case 'm':
        return NumTeammatesInCone(wid_dist_ratio,
                                  PlayerAbsolutePosition(side, num) + center, vert);
    case 't':
        return NumOpponentsInCone(wid_dist_ratio,
                                  PlayerAbsolutePosition(side, num) + center, vert);
    case 'b':
        return NumPlayersInCone(wid_dist_ratio,
                                PlayerAbsolutePosition(side, num) + center, vert);
    default:
        my_error("Bad which to NumPlayersInConeToPlayer");
        return 0;
    }
}

/*********************************************************************************/

Unum PositionInfo::ClosestTeammateTo(Vector p, Bool include_me)
{
    Unum ClosestPlayer = Unum_Unknown;
    float dist2, ClosestDist2 = SP_pitch_length * 2 * SP_pitch_length * 2;
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (!include_me && i == MyNumber)
            continue;
        if (TeammatePositionValid(i) && (dist2 = TeammateAbsolutePosition(i).dist2(p)) < ClosestDist2)
        {
            ClosestDist2 = dist2;
            ClosestPlayer = i;
        }
    }
    return ClosestPlayer;
}

/*********************************************************************************/

Unum PositionInfo::ClosestOpponentTo(Vector p)
{
    Unum ClosestPlayer = Unum_Unknown;
    float dist2, ClosestDist2 = SP_pitch_length * 2 * SP_pitch_length * 2;
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (OpponentPositionValid(i) && (dist2 = OpponentAbsolutePosition(i).dist2(p)) < ClosestDist2)
        {
            ClosestDist2 = dist2;
            ClosestPlayer = i;
        }
    }
    return ClosestPlayer;
}

/*********************************************************************************/

Vector PositionInfo::ClosestTeamlessPlayerPosition()
{
    if (num_teamless_players < 1)
        my_error("no teamless players");
    float tmp, closest_dist = HUGE;
    Vector position;
    for (int i = 0; i < num_teamless_players; i++)
    {
        if (!TeamlessPlayer[i]->pos_valid())
            my_error("teammless player not valid");
        if ((tmp = TeamlessPlayer[i]->get_dist()) < closest_dist)
        {
            position = TeamlessPlayer[i]->get_abs_pos();
            closest_dist = tmp;
        }
    }
    return position;
}

/*********************************************************************************/

Unum PositionInfo::ClosestTeammateToBall(Bool include_me)
{
    if (!BallPositionValid())
        my_error("can't do closest if ball unknown");
    return ClosestTeammateTo(BallAbsolutePosition(), include_me);
}

/*********************************************************************************/

Unum PositionInfo::ClosestOpponentToBall()
{
    if (!BallPositionValid())
        my_error("can't do closest if ball unknown");
    return ClosestOpponentTo(BallAbsolutePosition());
}

/*********************************************************************************/

float PositionInfo::ClosestTeammateToBallDistance(Bool include_me)
{
    Unum teammate = ClosestTeammateToBall(include_me);
    return BallAbsolutePosition().dist(TeammateAbsolutePosition(teammate));
}

/*********************************************************************************/

float PositionInfo::PlayerDistanceTo(char s, Unum n, Vector p)
{
    if (!PlayerPositionValid(s, n))
        my_error("can't get distance from invalid player");
    return PlayerAbsolutePosition(s, n).dist(p);
}

/*********************************************************************************/

float PositionInfo::PlayerDistance2To(char s, Unum n, Vector p)
{
    if (!PlayerPositionValid(s, n))
        my_error("can't get distance from invalid player");
    return PlayerAbsolutePosition(s, n).dist2(p);
}

/*********************************************************************************/

float PositionInfo::PlayerDistanceToLine(char s, Unum n, Line l)
{
    if (!PlayerPositionValid(s, n))
        my_error("can't get line distance from invalid player");
    return l.dist(PlayerAbsolutePosition(s, n));
}

/*********************************************************************************/

float PositionInfo::PlayerDistance2ToLine(char s, Unum n, Line l)
{
    if (!PlayerPositionValid(s, n))
        my_error("can't get line distance from invalid player");
    return l.dist2(PlayerAbsolutePosition(s, n));
}

/*********************************************************************************/

Unum PositionInfo::FurthestBackTeammate(Bool IncludeUs, Bool IncludeGoalie)
{
    Unum first = Unum_Unknown;
    float first_x = HUGE;
    Pnum tmp;
    for (Unum num = 1; num <= SP_team_size; num++)
    {
        if (!TeammatePositionValid(num))
            continue;
        if (!IncludeUs && num == MyNumber)
            continue;
        if (!IncludeGoalie && num == FP_goalie_number)
            continue;
        if (TeammateX(num) < first_x)
        {
            first = num;
            first_x = TeammateX(num);
        }
    }
    return first;
}

/*********************************************************************************/

Unum PositionInfo::FurthestBackOpponent()
{
    Unum first = Unum_Unknown;
    float first_x = HUGE;
    for (Unum num = 1; num <= SP_team_size; num++)
    {
        if (!OpponentPositionValid(num))
            continue;
        if (OpponentX(num) < first_x)
        {
            first = num;
            first_x = OpponentX(num);
        }
    }
    return first;
}

/*********************************************************************************/

Vector PositionInfo::PositionOfFurthestBackPlayer(Bool IncludeUs)
{
    Unum team = FurthestBackTeammate(IncludeUs);
    Unum opp = FurthestBackOpponent();
    if (team == Unum_Unknown && opp == Unum_Unknown)
        return Vector(0, 0); // no players found
    if (team == Unum_Unknown)
        return OpponentAbsolutePosition(opp);
    if (opp == Unum_Unknown)
        return TeammateAbsolutePosition(team);
    if (TeammateX(team) < OpponentX(opp))
        return TeammateAbsolutePosition(team);
    else
        return OpponentAbsolutePosition(opp);
}
/*********************************************************************************/

Unum PositionInfo::FurthestForwardTeammate(Bool IncludeUs)
{
    Unum first = Unum_Unknown;
    float first_x = -SP_pitch_length;
    for (Unum num = 1; num <= SP_team_size; num++)
    {
        if (!TeammatePositionValid(num))
            continue;
        if (!IncludeUs && num == MyNumber)
            continue;
        if (TeammateX(num) > first_x)
        {
            first = num;
            first_x = TeammateX(num);
        }
    }
    return first;
}

/*********************************************************************************/

Unum PositionInfo::FurthestForwardOpponent(Bool IncludeGoalie)
{
    Unum first = Unum_Unknown;
    Unum second = Unum_Unknown;
    float first_x = -SP_pitch_length;
    float second_x = -SP_pitch_length;
    for (Unum num = 1; num <= SP_team_size; num++)
    {
        if (!OpponentPositionValid(num))
            continue;
        if (OpponentX(num) > first_x)
        {
            second = first;
            second_x = first_x;
            first = num;
            first_x = OpponentX(num);
        }
        else if (OpponentX(num) > second_x)
        {
            second = num;
            second_x = OpponentX(num);
        }
    }

    if (!IncludeGoalie && first != Unum_Unknown &&
        TheirPenaltyArea.IsWithin(OpponentAbsolutePosition(first)))
        return second;

    return first;
}

/*********************************************************************************/

Vector PositionInfo::PositionOfFurthestForwardPlayer(Bool IncludeUs)
{
    Unum team = FurthestForwardTeammate(IncludeUs);
    Unum opp = FurthestForwardOpponent();
    if (team == Unum_Unknown && opp == Unum_Unknown)
        return Vector(0, 0); // no players found
    if (team == Unum_Unknown)
        return OpponentAbsolutePosition(opp);
    if (opp == Unum_Unknown)
        return TeammateAbsolutePosition(team);
    if (TeammateX(team) > OpponentX(opp))
        return TeammateAbsolutePosition(team);
    else
        return OpponentAbsolutePosition(opp);
}

/*********************************************************************************/

float PositionInfo::AngleBetweenClosestTwoOpponents(Vector p)
{
    Unum play1 = ClosestOpponentTo(p);
    Unum play2;
    float dist2, ClosestDist2 = Sqr(SP_pitch_length * 2);
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (i == play1)
            continue;
        if (OpponentPositionValid(i) &&
            (dist2 = OpponentAbsolutePosition(i).dist2(p)) < ClosestDist2)
        {
            ClosestDist2 = dist2;
            play2 = i;
        }
    }
    return fabs(GetNormalizeAngleDeg(OpponentAngleFromBody(play1) - OpponentAngleFromBody(play2)));
}

/*********************************************************************************/

Bool PositionInfo::InOwnPenaltyArea()
{
    if (MyConf() && InOwnPenaltyArea(MyPos()))
        return TRUE;
    else
        return FALSE;
}

/*********************************************************************************/

Bool PositionInfo::BallInOwnPenaltyArea()
{
    if (BallPositionValid() && InOwnPenaltyArea(BallAbsolutePosition()))
        return TRUE;
    else
        return FALSE;
}

/*********************************************************************************/

Bool PositionInfo::InOwnPenaltyArea(Vector p)
{
    if (p.x < MarkerX(RM_My_PC_Flag) && p.x > MarkerX(RM_My_Goal) &&
        fabs(p.y) < SP_penalty_area_width / 2.0)
        return TRUE;
    else
        return FALSE;
}

/*********************************************************************************/

Bool PositionInfo::FacingBackInOwnPA()
{
    /* neck facing back => can't see what's coming */
    if (InOwnPenaltyArea() && fabs(MyNeckGlobalAng()) > 90)
        return TRUE;
    else
        return FALSE;
}

/*********************************************************************************/

Bool PositionInfo::FacingBackNearOwnGoal()
{
    /* neck facing back => can't see what's coming */
    if (MyConf() && MarkerDistance(RM_My_Goal) < 25 && fabs(MyNeckGlobalAng()) > 90)
        return TRUE;
    else
        return FALSE;
}

/*********************************************************************************/

Bool PositionInfo::IsPointInBounds(float x, float y, float buffer)
{
    return (x > -SP_pitch_length / 2 + buffer && x < SP_pitch_length / 2 - buffer &&
            y > -SP_pitch_width / 2 + buffer && y < SP_pitch_width / 2 - buffer)
               ? TRUE
               : FALSE;
}

/*********************************************************************************/

Rectangle PositionInfo::ShiftRectangleInBounds(Rectangle *rect)
{
    float w = rect->Width();
    float h = rect->Height();

    Vector new_center = rect->Center();

    if (w > SP_pitch_length)
        my_error("rectangle too wide for field");
    if (h > SP_pitch_width)
        my_error("rectangle too high for field");

    if (rect->RightX() > SP_pitch_length / 2)
        new_center.x = SP_pitch_length / 2 - w / 2;
    if (rect->LeftX() < -SP_pitch_length / 2)
        new_center.x = -SP_pitch_length / 2 + w / 2;
    if (rect->BottomY() > SP_pitch_width / 2)
        new_center.y = SP_pitch_width / 2 - h / 2;
    if (rect->TopY() < -SP_pitch_width / 2)
        new_center.y = -SP_pitch_width / 2 + h / 2;

    return Rectangle(new_center, Vector(w, h));
}

/*********************************************************************************/

Vector PositionInfo::PositionToKickoffPosition(const Vector pos)
{

    Vector ko_pos = pos;

    if (ko_pos.x > 0 && fabs(ko_pos.y) < 9)
    {
        ko_pos.x = -1; /* Don't put the center forward behind the circle */
        ko_pos.y = 10;
    }

    if (fabs(ko_pos.x) <= 10 && ko_pos.y == 0)
    {
        if (KickOffMode == KO_Theirs)
        {
            ko_pos.x = -10;
            ko_pos.y = 0;
        }
        else
        {
            /* Put the center midfielder right at the ball */
            ko_pos.x = -CP_hardest_kick_ball_dist * Sin(45);
            // ko_pos.y = int_random(2) ? -(SP_kickable_area - 1) : SP_kickable_area-1;
            ko_pos.y = -CP_hardest_kick_ball_dist * Cos(45);
        }
    }
    else if (ko_pos.x > -1)
    {
        if (ko_pos.y > 0) /* Stagger the players     */
            ko_pos.y += ko_pos.x / 5;
        else
            ko_pos.y -= ko_pos.x / 5;
        ko_pos.x = -3; /* Always start on my side */
    }

    if (KickOffMode == KO_Theirs)
    { /* Stay out of the center circle */
        if (sqrt(ko_pos.x * ko_pos.x + ko_pos.y * ko_pos.y) < SP_free_kick_buffer)
            ko_pos.x = -10;
    }

    if (ko_pos.x > SP_pitch_length / 2)
        ko_pos.x = SP_pitch_length / 2;
    if (ko_pos.x < -SP_pitch_length / 2)
        ko_pos.x = -SP_pitch_length / 2;
    if (ko_pos.y > SP_pitch_width / 2)
        ko_pos.y = SP_pitch_width / 2;
    if (ko_pos.y < -SP_pitch_width / 2)
        ko_pos.y = -SP_pitch_width / 2;

    return ko_pos;
}

/*********************************************************************************/

/* consider_me : should I take myself into consideration when computing congestion */
float PositionInfo::Congestion(Vector pos, Bool consider_me)
{
    float congestion = 0;
    if (consider_me == TRUE && pos != MyPos())
        congestion = 1 / MyPos().dist2(pos);
    Vector player_pos;
    for (int i = 0; i < num_players; i++)
        /* Don't want to count a player in its own congestion measure */
        if (Player[i].pos_valid() && (player_pos = Player[i].get_abs_pos()) != pos)
            congestion += 1 / player_pos.dist2(pos);

    return congestion;
}

/*********************************************************************************/

float PositionInfo::TeammateCongestion(Unum teammate, Bool consider_me)
{
    if (!TeammatePositionValid(teammate))
        my_error("unknown teammate congestion");
    return Congestion(TeammateAbsolutePosition(teammate), consider_me);
}

/*********************************************************************************/

Unum PositionInfo::LeastCongestedTeammate()
{
    float least_congestion = Congestion(MyPos());
    Unum least_congested = MyNumber;
    Unum num;
    float temp;

    for (int i = 0; i < num_my_players; i++)
    {
        num = MyPlayer[i]->unum;
        if ((temp = Congestion(MyPlayer[i]->get_abs_pos(), TRUE)) < least_congestion)
        {
            least_congested = num;
            least_congestion = temp;
        }
    }
    return least_congested;
}

/*********************************************************************************/

Vector PositionInfo::LeastCongestedValidPointInRectangle(Rectangle *rect, Bool attract, Vector attract_point)
{
    int x_granularity = 5;
    int y_granularity = 5;

    float x_mesh = rect->Width() / (x_granularity + 1);
    float y_mesh = rect->Height() / (y_granularity + 1);

    float start_x = rect->LeftX() + x_mesh / 2;
    float start_y = rect->TopY() + y_mesh / 2;

    float x = start_x, y = start_y;

    float best_congestion = 1000;
    Vector best_point, point;
    float tmp;

    for (int i = 0; i < x_granularity; i++)
    {
        for (int j = 0; j < y_granularity; j++)
        {
            tmp = Congestion(point = Vector(x, y));
            if (attract)
                tmp -= 1 / point.dist(attract_point);
            if (tmp < best_congestion &&
                !OffsidePosition(point, MySide) &&
                IsPointInBounds(point, 5))
            {
                best_congestion = tmp;
                best_point = point;
            }
            y += y_mesh;
        }
        x += x_mesh;
        y = start_y;
    }

    if (best_congestion == 1000)
    {
        // my_error("No valid point in rectangle -- taking center %.1f %.1f",rect->LeftX(),rect->RightX());
        /* take the point out of the rectangle -- meaning no point was valid */
        best_point = rect->Center() + Vector(rect->Width(), 0);
    }

    return best_point;
}

/*********************************************************************************/

Vector PositionInfo::LeastCongestedValidPointForPassFromInRectangle(Rectangle *rect, Vector from, Bool attract, Vector attract_point)
{
    int x_granularity = 5;
    int y_granularity = 5;

    float x_mesh = rect->Width() / (x_granularity + 1);
    float y_mesh = rect->Height() / (y_granularity + 1);

    float start_x = rect->LeftX() + x_mesh / 2;
    float start_y = rect->TopY() + y_mesh / 2;

    float x = start_x, y = start_y;

    float best_congestion = 1000;
    Vector best_point, point;
    float tmp;

    for (int i = 0; i < x_granularity; i++)
    {
        for (int j = 0; j < y_granularity; j++)
        {
            tmp = Congestion(point = Vector(x, y));
            if (attract)
                tmp -= 1 / point.dist(attract_point);

            if (tmp < best_congestion &&
                !OffsidePosition(point, MySide) &&
                !NumOpponentsInCone(.6, point, from) && // was '&' in Paris
                IsPointInBounds(point, 5))
            {
                best_congestion = tmp;
                best_point = point;
            }
            y += y_mesh;
        }
        x += x_mesh;
        y = start_y;
    }

    if (best_congestion == 1000)
    {
        // my_error("No valid point in rectangle -- taking center %.1f %.1f",rect->TopY(),rect->BottomY());
        /* take the point out of the rectangle -- meaning no point was valid */
        best_point = rect->Center() + Vector(rect->Width(), 0);
    }

    return best_point;
}

/* -*- Mode: C++ -*- */

/* MemAction.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

#ifdef DEBUG_OUTPUT
#define DebugClear(x)
#else
#define DebugClear(x)
#endif

void ActionInfo::Initialize()
{

    Stored_Fastest_Teammate_Time = 0;
    Stored_Fastest_Opponent_Time = 0;

    for (int i = 1; i <= SP_team_size; i++)
    {
        TeamIntInfo[i] = new PlayerInterceptInfo;
        TeamIntInfo[i]->time = -1;
        OppIntInfo[i] = new PlayerInterceptInfo;
        OppIntInfo[i]->time = -1;
    }

    kick_in_progress = FALSE;

    InterceptLookahead = LA_Default;
    IntMinCyc = -1;
    IntMinCycTime = -1;

    HKTime = -1;
    HKStep = -1;
    HKStepNext = -1;
    HKrot = TURN_CW;
}

/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/

/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/

#ifdef DEBUG_OUTPUT
#define DebugInt(x)
#else
#define DebugInt(x)
#endif

/* only used for this player */
PlayerInterceptInfo
ActionInfo::CloseBallInterception(float max_pow, int max_lookahead,
                                  Vector vBallPos, Vector vBallVel)
{
    Vector vNewPos;
    PlayerInterceptInfo info;
    float dash_dist = max_pow * SP_dash_power_rate;
    info.dash_pow = max_pow;
    info.lookahead = max_lookahead;
    info.res = BI_None;

    vBallPos += vBallVel;
    vBallVel *= SP_ball_decay;
    /* we need to figure out the right dash power so that the ball ends up right in
       front of us. Like many things, this ends up as a LineCircleIntersect problem */
    Vector vMyPred = MyPredictedPosition();
    Ray rDash(vMyPred, MyBodyAng());
    Line lDash = LineFromRay(rDash);
    Vector sol1, sol2;

    int num_sol = LineCircleIntersect(lDash, SP_player_size + SP_ball_size + CP_collision_buffer,
                                      vBallPos, &sol1, &sol2);
    if (num_sol >= 1)
    {
        /* we'll make sure that the right answer is in sol1 */
        if (num_sol == 2)
        {
            /* we have to pick which point is better */
            if (fabs(GetNormalizeAngleDeg((vBallPos - sol2).dir() - MyBodyAng())) < 90)
            {
                sol1 = sol2; // sol2 is the right solution, put it in sol1
            }
            else if (!(fabs(GetNormalizeAngleDeg((vBallPos - sol1).dir() - MyBodyAng())) < 90))
            {
                my_error("CloseBallInterception: 1 ahead, neither solution looks good %.2f %.2f",
                         GetNormalizeAngleDeg((vBallPos - sol2).dir() - MyBodyAng()),
                         GetNormalizeAngleDeg((vBallPos - sol1).dir() - MyBodyAng()));
            }
        }

        /* now figure out the dash power based on that point */
        float dash_pow = vMyPred.dist(sol1) / SP_dash_power_rate;
        dash_pow = MinMax(SP_min_power, dash_pow, SP_max_power);
        if (!rDash.InRightDir(sol1))
            dash_pow = -dash_pow;

        if (vBallPos.dist(MyPredictedPosition(1, dash_pow)) < SP_kickable_area)
        {
            /* this works! */
            info.res = BI_CanChase;
            info.numCyc = 1;
            info.dash_pow_to_use = dash_pow;
            // this'll make go_to_point dash
            info.pos = vMyPred + Polar2Vector(signf(dash_pow) * dash_dist, MyBodyAng());
            LogAction5(70, "CloseBallInterception: One dash and we're there: %.2f to (%.2f, %.2f)",
                       dash_pow, sol1.x, sol1.y);
            return info;
        }
    }

    vBallPos += vBallVel;
    vBallVel *= SP_ball_decay;
    // now look two cycles ahead
    // try turn then dash
    float targ_ang = (vBallPos - MyPredictedPosition(2)).dir();
    if (fabs(targ_ang - MyBodyAng()) > CP_max_go_to_point_angle_err)
    {
        vNewPos = MyPredictedPositionWithTurn(targ_ang - MyBodyAng(), 2, max_pow);
        if (vNewPos.dist(vBallPos) < SP_kickable_area)
        {
            info.res = BI_CanChase;
            info.numCyc = 2;
            info.dash_pow_to_use = max_pow;
            // this'll make go_to_point turn
            // info.pos = MyPos() + Polar2Vector(dash_dist, targ_ang);
            info.pos = MyPredictedPosition() + Polar2Vector(dash_dist, targ_ang);
            LogAction2(70, "CloseBallInterception: Turn then dash and we're there");
            return info;
        }
    }
    // try two dashes
    vNewPos = MyPredictedPosition(2, max_pow);
    if (vNewPos.dist(vBallPos) < SP_kickable_area)
    {
        info.res = BI_CanChase;
        info.numCyc = 2;
        info.dash_pow_to_use = max_pow;
        // this'll make go_to_point dash
        // info.pos = MyPos() + Polar2Vector(2*dash_dist, MyBodyAng());
        info.pos = MyPredictedPosition() + Polar2Vector(dash_dist, MyBodyAng());
        LogAction2(70, "CloseBallInterception: Two dashes and we're there");
        return info;
    }
    return info;
}

/*****************************************************************************************/

/* does not set time field */
/* cyc inc is initially CP_intercept_step, but when an answer is found, we
   then bring cyc back a step, and go up by ones */
PlayerInterceptInfo
ActionInfo::ActiveCanGetThere(float max_pow, int max_lookahead,
                              Vector vBallPos, Vector vBallVel,
                              char side, Unum num,
                              Vector vPlayerPos, Vector vPlayerVel,
                              float fPlayerAng, int PlayerAngValid,
                              bool IsThisMe)
{
    float at_point_buffer = 1;
    PlayerInterceptInfo info;
    // Vector vPredPlayer = vPlayerPos + vPlayerVel;
    Vector vPredPlayer = vPlayerPos +
                         vPlayerVel * (SumInfGeomSeries(vPlayerVel.mod(), SP_player_decay));
    Vector vOldBallPos, vOldBallVel;
    float turn_angle;
    int cyc;
    int cyc_inc = (IsThisMe ? CP_my_intercept_step : CP_intercept_step);
    int max_cyc = (max_lookahead + cyc_inc - 1);
    max_cyc -= (max_cyc % cyc_inc);
    /* max_cyc is so that we don't miss an interception if CP_intercept_step is not
       1. For example, if CP_intercept_step is 5, max_look is 14, and we can
       intercept in 13, we would return no intercept if we just used max_lookahead */

    DebugInt(printf(" ACGT: BallPos.vel.mod: %f\n", vBallVel.mod()));

    info.dash_pow_to_use = max_pow;
    NormalizeAngleDeg(&fPlayerAng);

    /* we want to aim a little ahead of the ball, so advance it a little */
    for (int i = 0; i < CP_intercept_aim_ahead; i++)
    {
        vBallPos += vBallVel;
        vBallVel *= SP_ball_decay;
    }

    if (IsThisMe)
        LogAction4(140, "ActiveBallIntercept: %d %d", (int)max_pow, max_lookahead);

    for (cyc = 0; cyc <= max_cyc; cyc += cyc_inc)
    {

        if (!IsPointInBounds(vBallPos, -3))
        { /* expand the field by 3 meters so we don't give up to soon */
            DebugInt(printf("The ball will go out of bounds before we can get it\n"));
            break;
        }

        /* decide if we need to turn to ball */
        float ball_ang = (vBallPos - vPredPlayer).dir();
        Vector vEndSpot;
        /* SMURF - we should probably aim for ball 1 cycle ahead or something
           like that */
        // DebugInt(printf(" angle to exp ball pos: %f\n", AngleTo(vBallPos)));
        turn_angle = ball_ang - fPlayerAng;
        if (fabs(turn_angle) < CP_max_go_to_point_angle_err)
            turn_angle = 0.0;
        if (IsThisMe)
        {
            vEndSpot = MyPredictedPositionWithTurn(turn_angle, cyc, max_pow, (turn_angle != 0.0));
        }
        else
        {
            int run_cyc = cyc;
            if (PlayerAngValid)
            {
                if (turn_angle != 0.0)
                    run_cyc--;
                run_cyc = Max(0, run_cyc);
            }
            Vector PlayerDash =
                Polar2Vector(max_pow * SP_dash_power_rate, ball_ang);
            vEndSpot =
                PlayerPredictedPosition(side, num, run_cyc, PlayerDash);
        }

        float dist_to_ball_after = (vBallPos - vEndSpot).mod();
        /* if we can make it there */
        /* SMURF- is this too lenient? */
        if (dist_to_ball_after <= at_point_buffer ||
            (vEndSpot - vPredPlayer).mod() > (vBallPos - vPredPlayer).mod() + SP_kickable_area)
        {
            /* we can get to the ball! */
            /* OR we travelled far enough, but somehow missed the ball,
           return sucess */
            if (dist_to_ball_after <= at_point_buffer)
            {
                if (IsThisMe)
                    LogAction4(100, "Found a ball interception by being close (%.2f, %.2f)",
                               vBallPos.x, vBallPos.y);
                info.numCyc = cyc;
            }
            else
            {
                if (IsThisMe)
                    LogAction4(100, "Found a ball interception by going far (%.2f, %.2f)",
                               vBallPos.x, vBallPos.y);
                info.numCyc = cyc;
                // vBallPos += vBallVel; /* advance one spot for that turn*/
            }

            if (cyc_inc > 1 && cyc != 0)
            {
                /* we want the best answer- go back and go up by ones */
                if (IsThisMe)
                    LogAction2(100, "Found a ball interception, but goign back for accuracy");
                DebugInt(printf("Found answer, but going back for accuracy: %d\n", cyc));
                cyc -= cyc_inc;
                vBallPos = vOldBallPos;
                vBallVel = vOldBallVel;
                cyc_inc = 1;
                max_cyc = max_lookahead; // don;t need to go above this anymore
            }
            else
            {
                /* we want to try avoiding turning towards the ball for only a small savings
                   in time to intercept */
                if (IsThisMe && CP_no_turn_max_cyc_diff > -1 &&
                    turn_angle != 0.0 &&
                    (vBallVel.x >= FLOAT_EPS || vBallVel.y >= FLOAT_EPS))
                {
                    Ray rBall(vBallPos, vBallVel);
                    Ray rPlayer(vPredPlayer, fPlayerAng);
                    Vector int_pt;
                    if (rBall.intersection(rPlayer, &int_pt))
                    {
                        float dist = vEndSpot.dist(int_pt);
                        float num_cyc; /* the number of cycles extra it takes the ball to get to
                                this pos */
                        num_cyc = SolveForLengthGeomSeries(vBallVel.mod(), SP_ball_decay, dist);
                        LogAction3(90, "No turn interception: It takes %.2f extra cycles", num_cyc);
                        // if an answer less than 0 is given, the ball will never get there
                        if (num_cyc >= 0 &&
                            num_cyc <= CP_no_turn_max_cyc_diff)
                        {
                            /* use this target instead! */
                            LogAction4(70, "Using the new no turning interception point (%.2f, %.2f)",
                                       int_pt.x, int_pt.y);
                            info.res = BI_CanChase;
                            info.pos = int_pt;
                            return info;
                        } /* using no turn interseption */
                    }     /* there is an intersection */

                } /* no turn interseption */

                if (info.numCyc > max_lookahead)
                {
                    info.res = BI_Failure;
                }
                else
                {
                    info.res = BI_CanChase;
                    info.pos = vBallPos;
                }
                return info;
            }
        }

        /* update ball position estimate */
        vOldBallPos = vBallPos;
        vOldBallVel = vBallVel;
        for (int i = 0; i < cyc_inc; i++)
        {
            vBallPos += vBallVel;
            vBallVel *= SP_ball_decay;
        }

    } /* cycle loop */

    info.res = BI_Failure; // can't make it to ball before max_lookahead
    return info;
}

/*****************************************************************************************/

void ActionInfo::BallIntercept_active(float max_pow_to_use, int max_lookahead,
                                      char PlayerSide, Unum PlayerNum,
                                      PlayerInterceptInfo *pInfo)
{
    Vector PlayerPos;
    Vector PlayerVel;
    float PlayerAng;
    int AngValid = FALSE;
    Vector BallVel;

    pInfo->res = BI_None;

    if (!BallPositionValid())
    {
        my_error("BallIntercept_active: Can't get to ball if I don;t know where it is");
        pInfo->res = BI_Invalid;
        return;
    }

    if (!PlayerPositionValid(PlayerSide, PlayerNum))
    {
        my_error("BallIntercept_active: Can't give an answer if I don't know where player is");
        pInfo->res = BI_Invalid;
        return;
    }
    PlayerPos = PlayerAbsolutePosition(PlayerSide, PlayerNum);
    // DebugInt(cout << "PlayerPos: " << PlayerPos << endl);

    if (PlayerVelocityValid(PlayerSide, PlayerNum))
    {
        PlayerVel = PlayerAbsoluteVelocity(PlayerSide, PlayerNum);
    }
    else
    {
        PlayerVel = Vector(0, 0);
    }

    if (PlayerBodyAngleValid(PlayerSide, PlayerNum))
    {
        AngValid = TRUE;
        PlayerAng = PlayerAbsoluteBodyAngle(PlayerSide, PlayerNum);
    }
    else
        PlayerAng = 0;

    if ((PlayerPos - BallAbsolutePosition()).mod() <
        SP_kickable_area)
    {
        pInfo->res = BI_ReadyToKick;
        pInfo->numCyc = 0;
        pInfo->pos = PlayerPos;
        return;
    }

    if (BallVelocityValid())
        BallVel = BallAbsoluteVelocity();
    else
        BallVel = Vector(0, 0);

    DebugInt(printf("At BallIntercept_active  max_pow: %f, max_look: %d\n",
                    max_pow_to_use, max_lookahead));

    if (PlayerSide == MySide && PlayerNum == MyNumber)
        *pInfo = CloseBallInterception(max_pow_to_use, max_lookahead,
                                       BallAbsolutePosition(), BallVel);

    if (pInfo->res == BI_None)
        *pInfo =
            ActiveCanGetThere(max_pow_to_use, max_lookahead,
                              BallAbsolutePosition(), BallVel,
                              PlayerSide, PlayerNum,
                              PlayerPos, PlayerVel, PlayerAng, AngValid,
                              (PlayerSide == MySide && PlayerNum == MyNumber));
    else
        ; //{ printf("%d:%d Used Close Ball intercept\n",MyNumber,CurrentTime.t);}
}

/*****************************************************************************************/

PlayerInterceptInfo *ActionInfo::GetPlayerIntInfo(char side, Unum num)
{
    if (side == MySide)
        return TeamIntInfo[num];
    else if (side == TheirSide)
        return OppIntInfo[num];
    else
        my_error("bad side passed to GetPlayerIntInfo");
    return NULL;
}

/*****************************************************************************************/

PlayerInterceptInfo *ActionInfo::VerifyIntInfo(char side, Unum num, float dash_pow)
{
    PlayerInterceptInfo *pInfo = GetPlayerIntInfo(side, num);
    if (pInfo == NULL)
    {
        my_error("Bad side or number passed to VerifyIntInfo");
        return NULL;
    }

    int lookahead;
    switch (InterceptLookahead)
    {
    case LA_Default:
        lookahead = CP_max_int_lookahead;
        break;
    case LA_BestSoFar:
        lookahead =
            (IntMinCycTime == CurrentTime) ? (IntMinCyc) : CP_max_int_lookahead;
        break;
    default:
        lookahead = InterceptLookahead;
        break;
        break;
    }
    if (pInfo->time != CurrentTime || fabs((pInfo->dash_pow - dash_pow)) > FLOAT_EPS ||
        (pInfo->lookahead < lookahead && !IsSuccessRes(pInfo->res)))
    {
        /* set the info struct */
        DebugInt(printf("%d %d Data not current. Calling interception code\n", MyNumber, num));

        if (pInfo->time == CurrentTime && (pInfo->dash_pow - dash_pow) <= FLOAT_EPS &&
            (side != MySide || num != MyNumber))
            my_error("Recomputing %c %d because lookahead got bigger; old: %d\tnew: %d",
                     (int)side, (int)num, (int)pInfo->lookahead, (int)lookahead);

        /* let's do a real quick estimate to see if the player can make it there
           if player dist to ball > max ball dist will travel + max_dist we'll
           travel, then there's no way to get there */
        if (!PlayerPositionValid(side, num))
        {
            my_error("VerifyIntInfo: Can't give an answer if I don't know where player is");
            pInfo->res = BI_Invalid;
            return pInfo;
        }
        DebugInt(printf("Lookahead: %d\n", lookahead));
        float ball_travel = SumGeomSeries((BallVelocityValid() ? BallSpeed() : 0),
                                          SP_ball_decay, lookahead);
        float player_travel = SP_player_speed_max * lookahead;
        float play_ball_dist = (PlayerAbsolutePosition(side, num) -
                                BallAbsolutePosition())
                                   .mod();
        if (play_ball_dist > player_travel + ball_travel)
        {
            pInfo->time = CurrentTime;
            pInfo->dash_pow = dash_pow;
            pInfo->dash_pow_to_use = dash_pow;
            pInfo->lookahead = lookahead;
            pInfo->res = BI_Failure;
            DebugInt(printf("Interception: %d, %d Took shortcut to decide failure\n", MyNumber, num));
        }
        else
        {
            DebugInt(printf("Interception: %d, %d About to do actual calculation\n", MyNumber, num));
            BallIntercept_active(dash_pow, lookahead, side, num, pInfo);
            if (IsSuccessRes(pInfo->res))
                SetIntMinCyc(pInfo->numCyc);
            pInfo->time = CurrentTime;
            pInfo->dash_pow = dash_pow;
            pInfo->lookahead = lookahead;
        }
    }
    else if (IsSuccessRes(pInfo->res))
        SetIntMinCyc(pInfo->numCyc);

    return pInfo;
}

/*****************************************************************************************/

InterceptRes ActionInfo::PlayerInterceptionResult(char side, Unum num,
                                                  float dash_pow)
{
    return (VerifyIntInfo(side, num, dash_pow))->res;
}

/*****************************************************************************************/

Bool ActionInfo::PlayerInterceptionAble(char side, Unum num, float dash_pow)
{
    return IsSuccessRes((VerifyIntInfo(side, num, dash_pow))->res) ? TRUE : FALSE;
}

/*****************************************************************************************/

int ActionInfo::PlayerInterceptionNumberCycles(char side, Unum num,
                                               float dash_pow)
{
    PlayerInterceptInfo *pInfo = VerifyIntInfo(side, num, dash_pow);
    if (!IsSuccessRes(pInfo->res))
        my_error("Trying to get number of cycles on invalid result: %c%d %d",
                 side, num, pInfo->res);
    return pInfo->numCyc;
}

/*****************************************************************************************/

Vector ActionInfo::PlayerInterceptionPoint(char side, Unum num,
                                           float dash_pow)
{
    PlayerInterceptInfo *pInfo = VerifyIntInfo(side, num, dash_pow);
    if (!IsSuccessRes(pInfo->res))
        my_error("Trying to get interception point on invalid result: %c%d %d",
                 side, num, pInfo->res);
    return pInfo->pos;
}

/*****************************************************************************************/

float ActionInfo::PlayerInterceptionDashPower(char side, Unum num, float dash_pow)
{
    PlayerInterceptInfo *pInfo = VerifyIntInfo(side, num, dash_pow);
    if (!IsSuccessRes(pInfo->res))
        my_error("Trying to get interception dash power on invalid result: %c%d %d",
                 side, num, pInfo->res);
    return pInfo->dash_pow_to_use;
}

/*****************************************************************************************/

int ActionInfo::GetInterceptionMinCyc()
{
    if (IntMinCycTime != CurrentTime)
        return -1;
    else
        return IntMinCyc;
}

/*****************************************************************************************/

void ActionInfo::SetIntMinCyc(int newval)
{
    if (IntMinCycTime != CurrentTime)
    {
        IntMinCycTime = CurrentTime;
        IntMinCyc = newval;
    }
    else if (IntMinCyc > newval)
        IntMinCyc = newval;
}

/*****************************************************************************************/

void ActionInfo::SetInterceptionLookahead(int newval)
{
    if (newval > 0 || newval == LA_Default || newval == LA_BestSoFar)
    {
        if (IntMinCycTime == CurrentTime)
            DebugInt(cout << "Changing lookahead mid way through computations. Could be bad" << endl);
        InterceptLookahead = newval;
    }
    else
    {
        my_error("Trying to set InterceptLookahead to an invlaid value");
    }
}

/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/
/* Passive interception stuff */

int ActionInfo::GetClosestPointToBallPath(Vector *pPt, float *pNumCycles,
                                          Vector PlayerPos, Vector BallPos,
                                          Vector BallVel)
{
    if (fabs(BallVel.x) < FLOAT_EPS && fabs(BallVel.y) < FLOAT_EPS)
    {
        *pPt = BallPos;
        *pNumCycles = 0;
        return 1;
    }

    Ray rBallPath(BallPos, BallVel);
    ;

    *pPt = rBallPath.GetClosestPoint(PlayerPos);

    /* adjust point for sidelines */
    Rectangle field(Vector(0, 0), Vector(SP_pitch_length, SP_pitch_width));
    *pPt = AdjustPtToRectOnLine(*pPt, field, LineFromRay(rBallPath));

    /* Now let's reason about how far off we will be if we favor not turning */
    Vector no_turn_pt;
    if (rBallPath.intersection(Ray(MyPos(), MyBodyAng()), &no_turn_pt))
    {
        if (no_turn_pt.dist(*pPt) < CP_no_turn_max_dist_diff)
        {
            LogAction6(110, "BPI: using no turn interception, old: (%.1f, %.1f) new: (%.1f, %.1f)",
                       pPt->x, pPt->y, no_turn_pt.x, no_turn_pt.y);
            *pPt = no_turn_pt;
        }
    }

    /* compute the number of cycles to get here */
    *pNumCycles = 0;

    /* now get the number of cycles */
    Vector traj = *pPt - BallPos;
    DebugInt(cout << "Pt: " << *pPt << "\tBallVel: " << BallVel
                  << "\tBallPos: " << BallPos << "\ttraj: " << traj << endl);
    /* first decide if the ball is actually coming towards us */
    if (signf(traj.x) != signf(BallVel.x) ||
        signf(traj.y) != signf(BallVel.y))
    {
        DebugInt(printf("  GCPTBP: Ball is goign wrong way for closest intercept!\n"));
        return 0;
    }

    float trajDist = traj.mod();
    float velMod = BallVel.mod();
    float temp = trajDist / velMod * (SP_ball_decay - 1) + 1;
    if (temp < 0.0)
    {
        /* ball will never make it to closest point */
        /* SMURF - shoudl adjust for actual closest!!!! */
        DebugInt(printf("GCPTBP: Ball will never make it to closest point, adjusting\n"));
        *pPt = BallPos + traj * SumInfGeomSeries(velMod, SP_ball_decay) / traj.mod();
        *pNumCycles = SP_half_time; // just a big number
        return 1;
    }
    else
        *pNumCycles = log(temp) / log(SP_ball_decay);

    return 1;
}

/*****************************************************************************************/

void ActionInfo::VerifyBPIInfo()
{
    if (BPItime == CurrentTime)
        return;

    BPItime = CurrentTime;

    Vector BallVel;

    if (!MyConf())
    {
        my_error("Can't intercept if I don't know where I am");
        BPIvalid = FALSE;
        return;
    }

    if (!BallPositionValid())
    {
        my_error("Can't get to ball path if I don't know where it is");
        BPIvalid = FALSE;
        return;
    }

    if (BallKickable())
    {
        BPIvalid = TRUE;
        BPIable = TRUE;
        BPIdist = 0;
        BPIpoint = MyPos();
        BPIballcyc = 0;
        return;
    }

    if (BallVelocityValid())
        BallVel = BallAbsoluteVelocity();
    else
    {
        BPIvalid = TRUE;
        BPIable = TRUE;
        BPIdist = BallDistance();
        BPIpoint = BallAbsolutePosition();
        BPIballcyc = 0;
        return;
    }

    DebugInt(printf("\nTime: %d\n", CurrentTime.t));
    DebugInt(printf("At BallIntercept_passive\n"));

    int passRet;
    passRet = GetClosestPointToBallPath(&BPIpoint, &BPIballcyc, MyPos(),
                                        BallAbsolutePosition(), BallVel);
    DebugInt(printf("Passive Method: ret: %d\tx: %f\ty:%f\tcyc: %f\n",
                    passRet, BPIpoint.x, BPIpoint.y, BPIballcyc));
    if (passRet)
    {
        BPIvalid = TRUE;
        BPIable = TRUE;
        BPIdist = (BPIpoint - MyPos()).mod();
    }
    else
    {
        BPIvalid = TRUE;
        BPIable = FALSE;
    }

    return;
}

/*****************************************************************************************/

Vector ActionInfo::BallPathInterceptPoint()
{
    VerifyBPIInfo();
    if (!BPIvalid)
        my_error("Calling BallPathInterceptionPoint when info not valid?");
    return BPIpoint;
}

/*****************************************************************************************/

Bool ActionInfo::BallPathInterceptAmIThere(float buffer)
{
    VerifyBPIInfo();
    if (!BPIvalid)
        my_error("Calling BallPathInterceptionAmIThere when info not valid");
    return (BPIable && (MyPos() - BPIpoint).mod() <= buffer) ? TRUE : FALSE;
}

/*****************************************************************************************/

float ActionInfo::BallPathInterceptDistance()
{
    VerifyBPIInfo();
    if (!BPIable)
        my_error("Calling BallPathInterceptionDistance when I can't get get there");
    return BPIdist;
}

/*****************************************************************************************/

int ActionInfo::BallPathInterceptCyclesForBall()
{
    VerifyBPIInfo();
    if (!BPIable)
        my_error("Calling BallPathInterceptionCyclesForBall when I can't get get there");
    return (int)ceil(BPIballcyc);
}

/*****************************************************************************************/

Bool ActionInfo::BallPathInterceptCanIGetThere(float max_pow)
{
    VerifyBPIInfo();
    if (!BPIable)
        return FALSE;

    AngleDeg targAng = AngleToFromBody(BPIpoint);
    Vector myEnd;
    if (fabs(GetNormalizeAngleDeg(MyBodyAng() - targAng)) >
        CP_max_go_to_point_angle_err)
    {
        myEnd = MyPredictedPosition((int)ceil(BPIballcyc), max_pow);
    }
    else
    {
        myEnd = MyPredictedPositionWithTurn(targAng - MyBodyAng(),
                                            (int)ceil(BPIballcyc), max_pow);
    }

    return ((myEnd - MyPos()).mod() >= (BPIpoint - MyPos()).mod()) ? TRUE : FALSE;
}

/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/

float ActionInfo::VelAtPt2VelAtFoot(Vector pt, float targ_vel_at_pt)
{
    if (targ_vel_at_pt < FLOAT_EPS)
    {
        return SolveForFirstTermInfGeomSeries(SP_ball_decay, (pt - MyPos()).mod());
    }
    else
    {
        float ball_steps =
            SolveForLengthGeomSeries(targ_vel_at_pt, 1 / SP_ball_decay,
                                     (pt - MyPos()).mod());
        return targ_vel_at_pt * pow(1 / SP_ball_decay, ball_steps);
    }
}

/*****************************************************************************************/

/* looks at closeest opponent or teamless player */
/* SMURF: the teamless part is a hack */
KickMode ActionInfo::BestKickModeAbs(AngleDeg abs_ang)
{
    Unum closest = ClosestOpponent();
    if (NumTeamlessPlayers() > 0)
    {
        Vector teamless_pos = ClosestTeamlessPlayerPosition();
        if (closest == Unum_Unknown ||
            DistanceTo(teamless_pos) < OpponentDistance(closest))
            closest = Unum_Teamless;
    }

    if (closest == Unum_Unknown)
        return KM_HardestKick;
    int cyc_to_steal = EstimatedCyclesToSteal(closest);
    float targ_ang = abs_ang + signf(GetNormalizeAngleDeg(BallAngleFromBody() - abs_ang)) *
                                   (90 + CP_hardest_kick_ball_ang);
    float ang_diff = GetNormalizeAngleDeg(BallAngleFromBody() - targ_ang);
    NormalizeAngleDeg(&ang_diff);
    if (cyc_to_steal > fabs(ang_diff) / CP_time_for_full_rotation + CP_cycles_to_kick)
        return KM_HardestKick;
    // if (OpponentWithBall() != Unum_Unknown)
    if (cyc_to_steal <= 1)
        return KM_QuickestRelease;
    if (cyc_to_steal < CP_cycles_to_kick)
        return KM_Quickly;
    if (cyc_to_steal < CP_cycles_to_kick + 1) // time for a dash in KM_Hard
        return KM_Moderate;
    return KM_Hard;
}

/*****************************************************************************************/

/* returns estimated cycles for opponent to get the ball into his kickable
   area */
/* can handle Unum_Teamless SMURF: it's kind of a hack though */
int ActionInfo::EstimatedCyclesToSteal(Unum opp, Vector ball_pos)
{
    if (!BallKickable())
        my_error("EstimatedCyclesToSteal: shouldn't use this if the ball is not kickable");

    if (BallKickableForOpponent(opp))
    {
        LogAction2(110, "EstimatedCyclesToSteal: already kickable for opponent");
        return 0;
    }

    Vector targ = ball_pos;
    Vector pos;
    int cyc;
    if (opp == Unum_Teamless)
    {
        if (NumTeamlessPlayers() < 1)
            my_error("EstimatedCyclesToSteal: can't estimate teamless if there aren't any");
        pos = ClosestTeamlessPlayerPosition();
        targ -= (ball_pos - pos).SetLength(SP_kickable_area);
        cyc = (int)ceil(targ.dist(pos) / SP_player_speed_max);
    }
    else
    {
        if (!OpponentPositionValid(opp))
            my_error("EstimateCyclesToSteal: can't estimate if I don;t know where opponent is");
        pos = OpponentAbsolutePosition(opp);
        targ -= (ball_pos - pos).SetLength(SP_kickable_area);
        cyc = OpponentPredictedCyclesToPoint(opp, targ);
    }

    /* now decide if the player will have to dodge */
    if (!pos.ApproxEqual(targ))
    {
        Line oppLine = LineFromTwoPoints(pos, targ);
        Vector dodge_pos = oppLine.ProjectPoint(MyPos());
        dodge_pos += (pos - dodge_pos).SetLength(2.0 * SP_player_size);
        float dodge_dist = oppLine.dist(MyPos());
        if (dodge_dist < 2.0 * SP_player_size &&
            oppLine.InBetween(dodge_pos, pos, targ))
        {
            /* need to take into account a dodge */
            cyc += 2; // have to turn twice
            if (dodge_dist > 2 * SP_player_size - SP_player_speed_max)
                cyc += 1; // one dash will dodge us
            else
                cyc += 2; // it takes two dashes to dodge us
        }
    }

    return cyc;
}

/*****************************************************************************************/

/* this is not an exact function becuase we don't have a perfect mapping of
   ball speed/position to kick power.
   Basically, this function returns whether the ball will be further back but still
   kickable after a dash */
Bool ActionInfo::WillDashHelpKick(Vector pt, float dash_pow)
{
    if (!BallWillBeKickable(1, dash_pow, CP_kickable_buffer))
    {
        LogAction2(130, "WillDashHelpKick: ball will not be kickable");
        return FALSE;
    }

    /* we're going to assume that a collision is bad.
       but depending on how the ball is actually moving that could be good */
    if (WillDashBeCollision(dash_pow, CP_collision_buffer))
    {
        LogAction2(130, "WillDashHelpKick: collision");
        return FALSE;
    }

    /* if we're not facing genrally in the direction we want to kick it,
       dashing will probably not help */
    if (fabs(AngleToFromBody(pt)) > CP_max_dash_help_kick_angle)
    {
        LogAction2(130, "WillDashHelpKick: not facing");
        return FALSE;
    }

    AngleDeg curr_ang = BallAngleFromBody() - AngleToFromBody(pt);
    NormalizeAngleDeg(&curr_ang);
    Vector my_pred_pos = MyPredictedPosition(1, dash_pow);
    AngleDeg pred_ang =
        (BallPredictedPosition() - my_pred_pos).dir() -
        (pt - my_pred_pos).dir();
    NormalizeAngleDeg(&pred_ang);

    LogAction4(130, "WillDashHelpKick: curr: %.1f  pred: %.1f", curr_ang, pred_ang);

    return (fabs(pred_ang) > fabs(curr_ang)) ? TRUE : FALSE;
}

/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/

Bool ActionInfo::KickInProgress()
{
    /* need to have kicked last cycle.  Updates kick_in_progress_time */
    if (kick_in_progress && kick_in_progress_time == LastActionOpTime)
    {
        kick_in_progress_time = CurrentTime;
        return TRUE;
    }
    return FALSE;
}

/*****************************************************************************************/

void ActionInfo::StartKick(AngleDeg target_angle, KickMode mode, float target_vel, TurnDir rot)
{
    kick_in_progress = TRUE;
    start_kick_time = kick_in_progress_time = CurrentTime;
    kick_in_progress_abs_angle = GetNormalizeAngleDeg(target_angle + MyBodyAng());
    kick_in_progress_mode = mode;
    kick_in_progress_target_vel = target_vel;
    kick_in_progress_rotation = rot;
}

/*****************************************************************************************/

void ActionInfo::StartShot(AngleDeg target_angle, KickMode mode, TurnDir rot)
{
    StartKick(target_angle, mode, 2 * SP_ball_speed_max, rot);
}

/*****************************************************************************************/

void ActionInfo::StartPass(Unum target, float target_vel_at_dest, TurnDir rot)
{
    if (target == Unum_Unknown || !TeammatePositionValid(target))
        my_error("can't start this pass");

    team_passer = MyNumber;
    team_receiver = target;
    team_pass_time = CurrentTime;

    float target_vel = VelAtPt2VelAtFoot(TeammateAbsolutePosition(target), target_vel_at_dest);
    StartKick(TeammateAngleFromBody(target), KM_Moderate, target_vel, rot);
}

/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/

/* No reasoning about players being tired yet:
   if so, need to add dash_pow to the interception calls */

/* These functions are very computationally intensive */

/* the stored value does not include the goalie */
Unum ActionInfo::FastestTeammateToBall()
{
    if (!BallPositionValid())
        my_error("Need to know ball position to know fastest to it\n");

    Unum closest = ClosestTeammateToBall();
    if (!BallMoving() && !TeammateTired(closest))
        return closest;

    if (CurrentTime == Stored_Fastest_Teammate_Time)
        return Stored_Fastest_Teammate;

    ResetInterceptionMinCyc();
    SetInterceptionLookahead(LA_BestSoFar);

    Unum FastestPlayer = Unum_Unknown;
    int cycles, min_cycles = CP_max_int_lookahead + 1;

    for (int i = 1; i <= SP_team_size; i++)
    {
        if (TeammatePositionValid(i) && TeammateInterceptionAble(i) == TRUE &&
            (cycles = TeammateInterceptionNumberCycles(i)) < min_cycles &&
            (i != FP_goalie_number || CP_goalie))
        {
            min_cycles = cycles;
            FastestPlayer = i;
        }
    }

    Stored_Fastest_Teammate = FastestPlayer;
    Stored_Fastest_Teammate_Time = CurrentTime;

    return Stored_Fastest_Teammate;
}

/*****************************************************************************************/

Unum ActionInfo::FastestOpponentToBall()
{
    if (!BallPositionValid())
        my_error("Need to know ball position to know fastest to it\n");

    if (!BallMoving())
        return ClosestOpponentToBall();

    if (CurrentTime == Stored_Fastest_Opponent_Time)
        return Stored_Fastest_Opponent;

    ResetInterceptionMinCyc();
    SetInterceptionLookahead(LA_BestSoFar);

    Unum FastestPlayer = Unum_Unknown;
    int cycles, min_cycles = CP_max_int_lookahead + 1;
    for (int i = 1; i <= SP_team_size; i++)
    {
        if (OpponentPositionValid(i) && OpponentInterceptionAble(i) == TRUE &&
            (cycles = OpponentInterceptionNumberCycles(i)) < min_cycles)
        {
            min_cycles = cycles;
            FastestPlayer = i;
        }
    }

    Stored_Fastest_Opponent = FastestPlayer;
    Stored_Fastest_Opponent_Time = CurrentTime;

    return FastestPlayer;
}

/*****************************************************************************************/
Unum ActionInfo::BallPossessor()
{

    if (!BallPositionValid())
    {
        // my_error("BallPossesor: ball position not valid");
        return Unum_Unknown;
    }

    Unum num_with_ball = PlayerWithBall();
    if (num_with_ball != Unum_Unknown)
        return num_with_ball;

    Unum fastestTeammate, fastestOpponent;

    if (BallMoving())
    {
        int teamCycles, oppCycles;
        fastestOpponent = FastestOpponentToBall();
        fastestTeammate = FastestTeammateToBall();

        if (fastestTeammate == Unum_Unknown ||
            fastestOpponent == Unum_Unknown)
            return (fastestTeammate == Unum_Unknown ? -fastestOpponent : fastestTeammate);

        teamCycles = TeammateInterceptionNumberCycles(fastestTeammate);
        oppCycles = OpponentInterceptionNumberCycles(fastestOpponent);

        if (teamCycles + CP_possessor_intercept_space < oppCycles)
            return fastestTeammate;
        else if (oppCycles + CP_possessor_intercept_space < teamCycles)
            return -fastestOpponent;
    }
    else
    {
        fastestTeammate = ClosestTeammateToBall();
        fastestOpponent = ClosestOpponentToBall();

        if (fastestTeammate == Unum_Unknown ||
            fastestOpponent == Unum_Unknown)
            return (fastestTeammate == Unum_Unknown ? -fastestOpponent : fastestTeammate);

        /* we'll just ignore facing angles because they probably aren't right anwyay */;
        if (TeammateAbsolutePosition(fastestTeammate).dist(BallAbsolutePosition()) <
            OpponentAbsolutePosition(fastestOpponent).dist(BallAbsolutePosition()))
            return fastestTeammate;
        else
            return -fastestOpponent;
    }

    return Unum_Unknown;
}

/*****************************************************************************************/

char ActionInfo::TeamInPossession()
{
    switch (PlayMode)
    {
    case PM_Play_On:
        break;
    case PM_My_Kick_In:
    case PM_My_Corner_Kick:
    case PM_My_Kick_Off:
    case PM_My_Free_Kick:
    case PM_My_Goalie_Free_Kick:
    case PM_My_Offside_Kick:
    case PM_My_Goal_Kick:
        return MySide;
    case PM_Their_Kick_In:
    case PM_Their_Corner_Kick:
    case PM_Their_Goal_Kick:
    case PM_Their_Kick_Off:
    case PM_Their_Offside_Kick:
    case PM_Their_Free_Kick:
    case PM_Their_Goalie_Free_Kick:
        return TheirSide;
    default:
        break;
    }

    Unum player = BallPossessor();
    if (player > 0)
        return MySide;
    else if (player < 0)
        return TheirSide;
    else
        return '?';
}

/* -*- Mode: C++ -*- */

/* Memory.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

void Memory::Initialize()
{
    PlayerInfo::Initialize();
    PositionInfo::Initialize();
    ActionInfo::Initialize();
}

/* -*- Mode: C++ -*- */

/* client.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

void send_initialize_message();
void to_ori_pos();
void parse_initialize_message(char *);
Bool wait_for_signals(sigset_t *);
sigset_t init_handler();
void sigio_handler();
void sigalrm_handler();
void send_action();
void resend_last_action();

/* Global variables -- don't want to reallocate buffers each time */

sigset_t sigiomask, sigalrmask;

// Memory  Global_Mem;
// Memory *const Mem = &Global_Mem;

char recvbuf[MAXMESG];
char sendbuf[MAXMESG];

char *GLOBAL_sense_body_message = "(sense_body)";

int alrsigs_since_iosig = 0;

/****************************************************************************************/

int main(int argc, char *argv[])
{
    Mem = new Memory();

    if (Mem == NULL)
    {
        my_error("couldn't allocate Mem");
        exit(0);
    }

    Mem->GetOptions(argc, argv);

    Socket sock = init_connection(Mem->SP_host, Mem->SP_port);

    Mem->sock = &sock;

    if (Mem->sock->socketfd == -1)
    {
        std::cerr << "Can't open connection for player" << std::endl;
        exit(-1);
    }

    send_initialize_message();

    if (wait_message(recvbuf, Mem->sock) == 0)
        my_error("wait_message failed");

    parse_initialize_message(recvbuf);

    to_ori_pos();

    Mem->Initialize();

    sigset_t sigfullmask = init_handler();

    std::cout << "mem initialized" << std::endl;

    while (Mem->ServerAlive == TRUE && wait_for_signals(&sigfullmask))
    ;

    if (Mem->sock->socketfd != -1)
        close_connection(Mem->sock);

    printf("Shutting down player %d\n", Mem->MyNumber);
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void to_ori_pos() {
    const int init_pos[][2] = {
	// init team position
	{-5, 20},
	{-10, 0},
	{-5, -20},
	{-20, -28},
	{-20, -14},
	{-20, 14},
	{-20, 28},
	{-38, 15},
	{-32, 0},
	{-38, -15},
	{-45, 0}
    };
    char command[128];
    sprintf(command, "(move %d %d)", init_pos[Mem->MyNumber - 1][0], init_pos[Mem->MyNumber - 1][1]);
    send_message(command,Mem->sock);
}

/****************************************************************************************/

/* Send initialize message */
void send_initialize_message()
{
    if (Mem->IP_reconnect)
        sprintf(sendbuf, "(reconnect %s %d)", Mem->MyTeamName, Mem->IP_reconnect);
    else if (Mem->CP_goalie == TRUE && Mem->SP_version >= 4.00)
    {
        sprintf(sendbuf, "(init %s (version %.2f) (goalie))", Mem->MyTeamName, Mem->SP_version);
    }
    else
        sprintf(sendbuf, "(init %s (version %.2f))", Mem->MyTeamName, Mem->SP_version);

    if (send_message(sendbuf, Mem->sock) == -1)
        abort();
}

/****************************************************************************************/

/* Parse initialize message */
void parse_initialize_message(char *recvbuf)
{
    char mode[100];
    if (!(strncmp(recvbuf, "(init", 4)))
    {
        /* It's an init msg */
        sscanf(recvbuf, "(init %c %d %[^)]", &Mem->MySide, &Mem->MyNumber, mode);
        Mem->ServerAlive = TRUE;
    }
    else if (!(strncmp(recvbuf, "(reconnect", 4)))
    {
        /* It's a reconnect msg */
        sscanf(recvbuf, "(reconnect %c %[^)]", &Mem->MySide, mode);
        Mem->MyNumber = Mem->IP_reconnect;
        printf("reconnecting to %d on side %c!\n", Mem->MyNumber, Mem->MySide);
        Mem->ServerAlive = TRUE;
    }
    else
    {
        my_error("Didn't get an init message: '%s'", recvbuf);
        Mem->ServerAlive = FALSE;
    }

    if (Mem->CP_goalie && Mem->FP_goalie_number != Mem->MyNumber)
        my_error("goalie number inconsistent with me being goalie");

    if (!Mem->CP_goalie && Mem->FP_goalie_number == Mem->MyNumber)
        my_error("I should be the goalie");

    if (mode[0] == 'b')
    { /* Before_kick_off */
        Mem->SetPlayMode(PM_Before_Kick_Off);
        if (Mem->MySide == 'l')
            Mem->KickOffMode = KO_Mine;
        else
            Mem->KickOffMode = KO_Theirs;
    }
    else /* Act as if the game's in progress */
        Mem->SetPlayMode(PM_Play_On);
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

/* set time interval between the sensor receiving and command sending */
inline void set_timer()
{
    struct itimerval itv;
    itv.it_interval.tv_sec = 0;
    itv.it_interval.tv_usec = Mem->TimerInterval * 1000;
    itv.it_value.tv_sec = 0;
    itv.it_value.tv_usec = Mem->TimerInterval * 1000;
    setitimer(ITIMER_REAL, &itv, NULL);
}

inline void set_timer(int usec)
{
    struct itimerval itv;
    itv.it_interval.tv_sec = 0;
    itv.it_interval.tv_usec = Mem->TimerInterval * 1000;
    itv.it_value.tv_sec = 0;
    itv.it_value.tv_usec = usec;
    setitimer(ITIMER_REAL, &itv, NULL);
}

/****************************************************************************************/

sigset_t init_handler()
{
    sigemptyset(&sigalrmask);
    sigaddset(&sigalrmask, SIGALRM);
    sigemptyset(&sigiomask);
    sigaddset(&sigiomask, SIGIO);

    struct sigaction sigact;
    sigact.sa_flags = 0;
    sigact.sa_mask = sigiomask;

#ifdef Solaris
    sigact.sa_handler = (void (*)(int))sigalrm_handler;
#else
    sigact.sa_handler = (void (*)(int))sigalrm_handler;
#endif

    sigaction(SIGALRM, &sigact, NULL);
    sigact.sa_mask = sigalrmask;

#ifdef Solaris
    sigact.sa_handler = (void (*)(int))sigio_handler;
#else
    sigact.sa_handler = (void (*)(int))sigio_handler;
#endif

    sigaction(SIGIO, &sigact, NULL);
    set_timer();
    sigprocmask(SIG_UNBLOCK, &sigiomask, NULL);
    sigprocmask(SIG_UNBLOCK, &sigalrmask, NULL);

    sigset_t sigsetmask;
    sigprocmask(SIG_BLOCK, NULL, &sigsetmask); /* Get's the currently unblocked signals */
    return sigsetmask;
}

/****************************************************************************************/

/* suspend the process until one of the signals comes through */
/* could check for situation to kill client, return FALSE     */
/* i.e. too many actions with no sensory input coming in      */
Bool wait_for_signals(sigset_t *mask)
{
    sigsuspend(mask);
    return TRUE;
}

/****************************************************************************************/

/* SIGIO handler: receive and parse messages from server */
void sigio_handler()
{
    sigprocmask(SIG_BLOCK, &sigalrmask, NULL);
    int counter = 0;

    Time StartTime = Mem->CurrentTime;

    while (receive_message(recvbuf, Mem->sock) == 1)
    {
        Parse(recvbuf);
        counter++;
    }

    if (Mem->CurrentTime - StartTime > 1 && StartTime.s == 0 && Mem->CurrentTime.s == 0)
        my_error("Received several steps at once -- missing action ops!!! (%d %d)",
                 StartTime.t, StartTime.s);

    sigprocmask(SIG_UNBLOCK, &sigalrmask, NULL);

    alrsigs_since_iosig = 0;

    // if (counter>1) printf("Got %d messages\n",counter);
}

/****************************************************************************************/

/* SIGALRM handler: extract and send first command in commandlist */
void sigalrm_handler()
{
    sigprocmask(SIG_BLOCK, &sigiomask, NULL);

    if (Mem->LastInterruptTime != Mem->CurrentTime)
    {
        if (!Mem->ClockStopped && Mem->CurrentTime - 1 != Mem->LastInterruptTime)
            my_error("Missed a cycle??");
        if (!Mem->ClockStopped && Mem->InterruptsThisCycle < Mem->CP_interrupts_per_cycle - 1)
            my_error("Only %d interrupts last cycle", Mem->InterruptsThisCycle);
        Mem->LastInterruptTime = Mem->CurrentTime;
        Mem->InterruptsThisCycle = 0;
        // cout << endl;
    }
    Mem->InterruptsThisCycle++;

    // cout << ".";

    /* Don't act until near the end of a cycle */
    /* there's some leeway in case there aren't enough interrupts in the cycle */
    if (!Mem->ClockStopped && Mem->CP_interrupts_per_cycle - Mem->InterruptsThisCycle >
                                  Mem->CP_interrupts_left_to_act)
        return;

    if (Mem->ClockStopped)
        Mem->StoppedClockMSec += Mem->TimerInterval;

    if (alrsigs_since_iosig++ > Mem->CP_interrupts_per_cycle * 20)
    {
        Mem->ServerAlive = FALSE;
        return;
    }

    /* If a sight is definitely coming every cycle, don't act until getting the sight */
    /* Don't wait if we're in transition to a longer sight interval                   */
    if (Mem->MySightInterval() < Mem->SP_simulator_step && Mem->LastSightTime < Mem->CurrentTime &&
        !((Mem->ChangeView.valid() || Mem->ChangeView.valid(Mem->CurrentTime - 1)) &&
          (Mem->ChangeView.width > Mem->ViewWidth || Mem->ChangeView.qual > Mem->ViewQuality)))
    {
        Mem->LogAction4(200, "Waiting for sight... (%d %d)",
                        Mem->ChangeView.valid(), Mem->ChangeView.valid(Mem->CurrentTime - 1));
        return;
    }

    if (Mem->CurrentTime > Mem->LastActionOpTime)
    {

        if (!Mem->ClockStopped && Mem->CurrentTime - 1 != Mem->LastActionOpTime && Mem->LastActionOpTime != 0)
            my_error("Missed a cycle!!  (%d %d)", Mem->LastActionOpTime.t, Mem->LastActionOpTime.s);

        if (Mem->NewSight)
            Mem->FirstActionOpSinceLastSight = TRUE;

        /*if ( 0 && Mem->MyCurrentFormationType() != FT_433 ){ my_error("?? pre-update %d\n"); dump_core("dump"); }*/

        Mem->update();
        behave();

        /*if ( 0 && Mem->MyCurrentFormationType() != FT_433 ){ my_error("?? post-behave %d\n"); dump_core("dump"); }*/

        Mem->LastActionOpTime = Mem->CurrentTime;
        Mem->FirstActionOpSinceLastSight = FALSE;
    }

    /*if ( 0 && Mem->MyCurrentFormationType() != FT_433 ){ my_error("?? pre-action %d\n"); dump_core("dump"); }*/

    /* Whether or not to wait between sending network packets is an interesting decsision.
       In versions of the server after 5.23, the server reads *all* commands off a socket at
       at every time step, so we could try to send all of our commands as soon as they are
       ready. However, on an ethernet, this can lead to lots of collisions and such, so it
       may degrade network performance
       To send everything without waiting, comment in this next line */
    //#define SEND_ALL_AT_ONCE

    /* the server now accepts multiple commands together  (after 5.23)
    if (0 && Mem->TooSoonForAnotherSend()) {
      Mem->LogAction2(200, "It's too soon to send another command. Waiting");
    } else {
    */

    if (Mem->Action->valid())
    {
        send_action();
    }
#ifndef SEND_ALL_AT_ONCE
    else
#endif
        if (Mem->ResendNeeded())
    {
        resend_last_action();
    }
#ifndef SEND_ALL_AT_ONCE
    else
#endif
        if (Mem->TurnNeck.valid())
    {
        if (Mem->TurnNeck.time < Mem->CurrentTime - 1)
            my_error("old turn_neck");

        send_message(Mem->TurnNeck.command, Mem->sock);
        Mem->turn_necks++;
        Mem->TurnNeck.type = CMD_none; /* so it's no longer valid */
    }
#ifndef SEND_ALL_AT_ONCE
    else
#endif
        if (Mem->ChangeView.valid())
    {
        if (Mem->ChangeView.time < Mem->CurrentTime - 1)
            my_error("old change_view");

        send_message(Mem->ChangeView.command, Mem->sock);
        Mem->ChangeView.type = CMD_none; /* so it's no longer valid */
    }
#ifndef SEND_ALL_AT_ONCE
    else
#endif
        if (Mem->SP_sense_body_step > Mem->SP_simulator_step)
    {
        /* only if we won't get a sense_body each cycle by default */
        my_error("Sending sense_body");
        send_message(GLOBAL_sense_body_message, Mem->sock);
    }

    sigprocmask(SIG_UNBLOCK, &sigiomask, NULL);
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

/* insert turn/dash/kick commands in commandlist */

void turn(AngleDeg ang)
{
    NormalizeAngleDeg(&ang);

    /* turn so that the actual turn is the desired turn */
    /* pos.rotate( ang/(1.0 + SP_inertia_moment * MySpeed()) ); */
    if (Mem->MyVelConf())
        ang *= (1 + Mem->SP_inertia_moment * Mem->MySpeed());

    if (ang > Mem->SP_max_moment)
        ang = Mem->SP_max_moment;
    if (ang < Mem->SP_min_moment)
        ang = Mem->SP_min_moment;

    if (ang < .1 && ang > -.1)
    {
        Mem->Action->type = CMD_none;
        return; /* No turn           */
    }

    Mem->Action->type = CMD_turn;
    Mem->Action->power = 0;
    Mem->Action->angle = ang;
    Mem->Action->time = Mem->CurrentTime;

    if (Mem->TurnNeck.time == Mem->CurrentTime)
    { /* Already did a turn_neck     */
        /* adjust as much as possible for the turn */
        Mem->TurnNeck.angle -= ang;
        if (Mem->MyNeckRelAng() + Mem->TurnNeck.angle < Mem->SP_min_neck_angle)
            Mem->TurnNeck.angle = Mem->SP_min_neck_angle - Mem->MyNeckRelAng();
        if (Mem->MyNeckRelAng() + Mem->TurnNeck.angle > Mem->SP_max_neck_angle)
            Mem->TurnNeck.angle = Mem->SP_max_neck_angle - Mem->MyNeckRelAng();
    }

    sprintf(Mem->Action->command, "(turn %.2f)", ang);
    Mem->LogAction3(150, "turn %f", ang);
}

/****************************************************************************************/

void dash(float power)
{
    if (Mem->PlayMode == PM_Before_Kick_Off)
        return;

    if (power > Mem->SP_max_power)
        my_error("Can't dash that fast: %.1f", power);
    if (power < Mem->SP_min_power)
        my_error("Can't dash that 'slow': %.1f", power);

    /* Factor for stamina--don't dash more than stamina or more than necessary to get you to max speed */
    Mem->VerifyDash(&power);

    if (fabs(power) < 1)
    {
        Mem->Action->type = CMD_none;
        return; /* No dash           */
    }

    Mem->Action->type = CMD_dash;
    Mem->Action->power = power;
    Mem->Action->angle = 0;
    Mem->Action->time = Mem->CurrentTime;

    sprintf(Mem->Action->command, "(dash %.2f)", power);
    Mem->LogAction3(150, "dash %f", power);
}

/****************************************************************************************/

void kick(float power, AngleDeg dir)
{
    if (!(Mem->BallKickable()))
        my_error("Can't kick a ball that's too far away");

    if (Mem->PlayMode == PM_Before_Kick_Off)
        return;

    if (power > Mem->SP_max_power)
        my_error("Can't kick that hard");
    if (power < 0)
        my_error("Can't kick < 0");
    NormalizeAngleDeg(&dir);

    Mem->Action->type = CMD_kick;
    Mem->Action->power = power;
    Mem->Action->angle = dir;
    Mem->Action->time = Mem->CurrentTime;

    sprintf(Mem->Action->command, "(kick %.2f %.2f)", power, dir);
    Mem->LogAction4(150, "kick %f %f", power, dir);
}

/****************************************************************************************/

void goalie_catch(AngleDeg dir)
{
    if (!(Mem->BallCatchable()))
        my_error("Can't catch a ball that's too far away");
    if (!Mem->CP_goalie)
        my_error("Only goalies can catch");

    if (Mem->PlayMode == PM_Before_Kick_Off)
        return;

    NormalizeAngleDeg(&dir);

    Mem->Action->type = CMD_catch;
    Mem->Action->power = 0;
    Mem->Action->angle = dir;
    Mem->Action->time = Mem->CurrentTime;

    sprintf(Mem->Action->command, "(catch %.2f)", dir);
    Mem->LogAction3(150, "catch %f", dir);
}

/****************************************************************************************/

void move(float x, float y)
{
    if (!(Mem->PlayMode == PM_Before_Kick_Off ||
          (Mem->CP_goalie && Mem->PlayMode == PM_My_Goalie_Free_Kick)))
        my_error("Can only move in before kickoff mode (or after goalie catch)");

    /* Perhaps here convert to a position on the field */
    if (fabs(y) > Mem->SP_pitch_width / 2 || x > 0 || x < -Mem->SP_pitch_length / 2)
        my_error("Must move to a place on the pitch");

    if (Mem->PlayMode == PM_My_Goalie_Free_Kick && !Mem->OwnPenaltyArea.IsWithin(Vector(x, y)))
        my_error("Must move to a place within penalty area");

    Mem->Action->type = CMD_move;
    Mem->Action->x = x;
    Mem->Action->y = y;
    Mem->Action->time = Mem->CurrentTime;

    sprintf(Mem->Action->command, "(move %.2f %.2f)", x, y);
    Mem->LogAction4(150, "move %f %f", x, y);
}

/****************************************************************************************/

void disconnect()
{
    Mem->Action->type = CMD_bye;
    Mem->Action->time = Mem->CurrentTime;

    sprintf(Mem->Action->command, "(bye)");
}

/****************************************************************************************/

void turn_neck(AngleDeg ang)
{
    NormalizeAngleDeg(&ang);

    if (ang == 0)
    {
        Mem->LogAction2(150, "Ignoring turn_neck 0");
        return;
    }

    if (ang > Mem->SP_max_neck_moment)
        ang = Mem->SP_max_neck_moment;
    if (ang < Mem->SP_min_neck_moment)
        ang = Mem->SP_min_neck_moment;

    if (Mem->MyNeckRelAng() + ang > Mem->SP_max_neck_angle)
    {
        ang = Mem->SP_max_neck_angle - Mem->MyNeckRelAng();
        my_error("Can't turn neck that much");
    }

    if (Mem->MyNeckRelAng() + ang < Mem->SP_min_neck_angle)
    {
        ang = Mem->SP_min_neck_angle - Mem->MyNeckRelAng();
        my_error("Can't turn neck that little");
    }

    Mem->TurnNeck.type = CMD_turn_neck;
    Mem->TurnNeck.power = 0;
    Mem->TurnNeck.angle = ang;
    Mem->TurnNeck.time = Mem->CurrentTime;

    sprintf(Mem->TurnNeck.command, "(turn_neck %.2f)", ang);
    Mem->LogAction3(150, "turn_neck %f", ang);
}

/****************************************************************************************/

void change_view(Vqual qual, Vwidth width)
{
    if (qual == Mem->ViewQuality && width == Mem->ViewWidth)
        return; /* my_error("Nothing to change about view"); */

    Mem->ChangeView.type = CMD_change_view;
    Mem->ChangeView.qual = qual;
    Mem->ChangeView.width = width;
    Mem->ChangeView.time = Mem->CurrentTime;

    char qual_string[10], width_string[10];
    switch (qual)
    {
    case VQ_High:
        sprintf(qual_string, "high");
        break;
    case VQ_Low:
        sprintf(qual_string, "low");
        break;
    }
    switch (width)
    {
    case VW_Narrow:
        sprintf(width_string, "narrow");
        break;
    case VW_Normal:
        sprintf(width_string, "normal");
        break;
    case VW_Wide:
        sprintf(width_string, "wide");
        break;
    }

    sprintf(Mem->ChangeView.command, "(change_view %s %s)", width_string, qual_string);
    Mem->LogAction4(150, "change_view %s %s", width_string, qual_string);
}

/****************************************************************************************/

void send_action()
{
    if (!(Mem->Action->valid(Mem->CurrentTime)))
        my_error("Old action %d %d", Mem->Action->time.t, Mem->Action->time.s);

    send_message(Mem->Action->command, Mem->sock);

    switch (Mem->Action->type)
    {
    case CMD_kick:
        Mem->kicks++;
        Mem->GetBall()->set_past_kick(Mem->LastActionPower(), Mem->LastActionAngle(),
                                      Mem->CurrentTime);
        break;
    case CMD_dash:
        Mem->dashes++;
        break;
    case CMD_turn:
        Mem->turns++;
        break;
    default:;
    }

    Command *tmp = Mem->LastAction;
    Mem->LastAction = Mem->Action;
    Mem->Action = tmp;

    Mem->Action->type = CMD_none; /* So it's invalid */
    Mem->NewAction = TRUE;
}

/****************************************************************************************/

void resend_last_action()
{
    if (Mem->LastActionType() == Mem->ResendType)
    {
        my_stamp;
        printf("resending\n");
        send_message(Mem->LastAction->command, Mem->sock);

        switch (Mem->LastActionType())
        {
        case CMD_kick:
            Mem->kicks++;
            break;
        case CMD_dash:
            Mem->dashes++;
            break;
        case CMD_turn:
            Mem->turns++;
            break;
        default:;
        }
    }
    else
        my_error("last action isn't a %d", Mem->ResendType);

    Mem->RequestResend = FALSE;
}

/* -*- Mode: C++ -*- */

/* kick.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

#ifdef DEBUG_OUTPUT
#define DebugKick(x)
#define DebugKick2(x)
#else
#define DebugKick(x)
#define DebugKick2(x)
#endif

/* these PatsTest_* functions are just that.
   They are temporary function that can be used to see some specific kicking
   behavior.
   What they do isn't documented anywhere but the source code below */

void PatsTest_static()
{
    DebugKick(printf("\nTime: %d\n", Mem->CurrentTime.t));
    if (!Mem->MyConf())
        scan_field_with_body();
    else if (!Mem->BallPositionValid())
        face_neck_to_ball();
    else
    {
        if (Mem->CurrentTime.t % 100 < 10 && Mem->BallKickable())
            smart_kick_hard_abs(45, KM_Moderate);
        else
        {
            Bool res = go_to_static_ball(45);
            std::cout << "result: " << res << std::endl;
        }
    }
}

void PatsTest_conv(void)
{
    Vector g, r, g2;
    for (float i = -10.0; i <= 10.0; i += 1.0)
    {
        for (float j = -10.0; j <= 10.0; j += 1.0)
        {
            g = Vector(i, j);
            r = g.Global2Relative(-3.0, 37);
            g2 = r.Relative2Global(-3.0, 37);
            if (fabs(g.x - g2.x) > FLOAT_EPS ||
                fabs(g.y - g2.y) > FLOAT_EPS)
                printf("Error i:%f\tj: %f\n", i, j);
        }
    }
    exit(1);
}

void PatsTest_turn(void)
{
    KickToRes res;
    static int ang = 180;
    static int wait_count = 0;

    /* SMURF - there is something wrong with Time class!
    if (Mem->CurrentTime - Mem->LastActionTime < Mem->CP_kick_time_space)
      return;
      */
    // TurnKickCommand com;
    if (Mem->CurrentTime.t % 100 < 20)
        return;

    res = KT_None;
    if (Mem->BallKickable())
    {

        if (wait_count <= 0)
        {
            DebugKick(printf("About to call turnball_kick\n"));
            res = TurnballTo((AngleDeg)ang - Mem->MyBodyAng(), TURN_CLOSEST);
        }
        else
        {
            DebugKick(printf("waiting...\n"));
            wait_count--;
        }
    }
    else
    {
    }

    if (res == KT_LostBall)
        return;
    else if (res == KT_Success)
    {
        DebugKick(printf("SUCESS-play_with_ball!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"));
        DebugKick(printf("ang: %d", ang));
        ang = (90 + ang) % 360;
        DebugKick(printf("new ang: %d\n", ang));
    }
    else if (res == KT_DidNothing)
        return;
    return;
}

void PatsTest_kick()
{
    if (!strcmp(Mem->MyTeamName, "CMUnited"))
    {
        // float targ_vel;
        Vector pass_targ;

        fflush(stdout);

        change_view(VW_Narrow);

        if (!Mem->BallKickable())
        {
            // scan_field();
            return;
        }

        /*    if((int)(Mem->CurrentTime.t) % 50 < 5 ) {
          return;
        }*/

        DebugKick(printf("\nTime: %d\n", Mem->CurrentTime.t));

        /*
        targ_vel = 1.0 + (Mem->CurrentTime.t / 50)*.1;
        printf("Target vel is: %f\n", targ_vel);*/
        pass_targ = Vector(-20, -20 + (Mem->CurrentTime.t / 50) * 5);
        // cout << "Pass Target: " << pass_targ << endl;

        // DebugKick(printf("the (simulated) ball velocity is: %g\n",)
        //	 Mem->BallAbsoluteVelocity().mod());
        /*
        DebugKick(printf("My angle: %f\n", Mem->MyAng()));
        if (fabs(Mem->MyAng() - 90) > 5) {
          DebugKick(printf("Turning to 90\n"));
          //turn(90-Mem->MyAng());
          turn(Mem->MarkerAngle(Mem->RM_RC_Flag));
          return;
        }
        */
        if (Mem->BallKickable())
        {
            /*if (step == 4)
          if (step == 3) {
          step = 0;
          turn(90 - Mem->MyAng());
          } 	*/
            smart_kick_hard_abs(180, KM_Moderate);
            // smart_kick_hard_abs(180, KM_HardestKick);
            // smart_kick_hard_abs(180, KM_Moderate, targ_vel);
            // smart_pass(pass_targ);
        }
    }
    else
    {
        static int FirstTime = TRUE;
        const int ydist = 2;
        //    if (Mem->MyNumber == 1)
        //      test_go_to_point(Vector(-.2, -1), .5, 50);
        if (!FirstTime)
            return;
        switch (Mem->MyNumber)
        {
        case 1:
            move(-.2, -1);
            break;
        case 2:
            move(-20, ydist);
            break;
        case 3:
            move(-30, -ydist);
            break;
        case 4:
            move(-40, ydist);
            break;
        case 5:
            move(-50, -ydist);
            break;
        }
        FirstTime = FALSE;
    }

    return;
}

void PatsTest_kick2()
{
    fflush(stdout);

    if (Mem->ViewWidth != VW_Wide)
        change_view(VW_Wide);

    static int WasBallKickable;
    DebugKick(printf("Time: %d\n", Mem->CurrentTime.t));
    // DebugKick(printf(" Ball Distance: %f\n", Mem->BallDistance() ));
    // DebugKick(cout << " MyPos: " << Mem->MyPos() << endl));
    if (!Mem->BallPositionValid() || !Mem->BallKickable())
    {
        WasBallKickable = FALSE;

        DebugKick(printf("chasing ball\n"));
        if (Mem->BallPositionValid())
        {
            if (fabs(Mem->BallAngleFromBody()) > Mem->CP_KickTo_err)
            {
                DebugKick(printf("turning to face ball\n"));
                DebugKick(printf("the angle is: %f\n", Mem->BallAngle()));
                turn(Mem->BallAngleFromBody());
            }
            else
            {
                DebugKick(printf("dashing\n"));
                float power = Mem->BallDistance() * 40;
                if (power > 40)
                    power = 40;
                dash(power);
            }
        }
        else
        {
            DebugKick(printf("turning randomly\n"));
            turn(60);
        }
        return;
    }

    // DebugKick(printf("the (simulated) ball velocity is: %g\n",)
    //	 Mem->BallAbsoluteVelocity().mod());

    if (Mem->BallKickable())
    {
        if (!smart_kick_hard_abs(0, KM_HardestKick))
            printf("test_kick_hard: UhOh, something bad happened\n");
        return;
    }

    return;
}

int DoTurnKickCommand(TurnKickCommand com)
{
    if (com.time != Mem->CurrentTime)
    {
        my_error("DoTurnKickCommand- told to do command not set this cycle");
        return 0;
    }

    switch (com.type)
    {
    case CMD_dash:
        DebugKick(printf("DoTurnKickCommand: dash\n"));
        dash(com.power);
        break;
    case CMD_turn:
        DebugKick(printf("DoTurnKickCommand: turn\n"));
        turn(com.angle);
        break;
    case CMD_kick:
        DebugKick(printf("DoTurnKickCommand: kick\n"));
        kick(com.power, com.angle);
        break;

    default:
        my_error("PatInfo::DoTurnKickCommand- unimplemented type!");
        return 0;
    }

    if (com.turn_neck)
    {
        turn_neck(com.turn_neck_angle);
    }
    return 1;
}

/* decides if we can kick staright to the decired point around the player
   without a collision.
   (EndDist, dir) is a relative polar vector for the ball's final position
   closeMarg is what the radius of the player is considered to be */
int is_straight_kick(float dir, float EndDist, float closeMarg)
{
    Vector btraj = Polar2Vector(EndDist, dir) - Mem->BallRelativeToBodyPosition();
    float ang;
    float dist;
    int res;

    DebugKick(printf("    isStriaght ball abs pos mod: %f\t dir: %f\n",
                     Mem->BallAbsolutePosition().mod(), Mem->BallAbsolutePosition().dir()));
    DebugKick(printf("    isStriaght ball rel pos mod: %f\t dir: %f\n",
                     Mem->BallRelativePosition().mod(), Mem->BallRelativePosition().dir()));
    DebugKick(printf("    isStriaght btraj mod: %f\t dir: %f\n", btraj.mod(), btraj.dir()));
    /* Apply the law of cosines to the anle formed by the player's center(A),
       the ball's current position(B), and the ball target position(C).
       The angle calculated is ABC */
    ang = ACos((Sqr(EndDist) - Sqr(btraj.mod()) -
                Sqr(Mem->BallDistance())) /
               (-2 * Mem->BallDistance() * btraj.mod()));
    DebugKick(printf("   isStraight ang: %f\n", ang));
    if (fabs(ang) > 90)
    {
        DebugKick(printf("    isStraight: Obtuse!\n"));
        Mem->LogAction2(120, "is_straight_kick: obtuse angle");
        return 1; /* obtuse angle implies definately straight */
    }
    /* get the height of the triangle, ie how close to the player the ball will
       go */
    dist = Sin(ang) * Mem->BallDistance();
    DebugKick(printf("   isStraight dist: %f\n", dist));
    Mem->LogAction3(120, "is_straight_kick: %f", dist);
    res = (fabs(dist) > closeMarg);
    return (res);
}

/* picks a rotation (CW or CCW) to turnball based on the nearest opponent or
   teamless player */
TurnDir RotToAvoidOpponent(float abs_dir)
{
    TurnDir rot;
    float dist = HUGE;
    AngleDeg opp_ang;
    Unum opp = Mem->ClosestOpponent();
    AngleDeg ball_ang = Mem->BallAngleFromBody();
    if (opp != Unum_Unknown)
    {
        dist = Mem->OpponentDistance(opp);
        opp_ang = Mem->OpponentAngleFromBody(opp);
    }
    if (Mem->NumTeamlessPlayers() > 0)
    {
        Vector pos = Mem->ClosestTeamlessPlayerPosition();
        float d = Mem->DistanceTo(pos);
        if (d < dist)
        {
            dist = d;
            opp_ang = Mem->AngleToFromBody(pos);
        }
    }
    if (dist < Mem->CP_turnball_opp_worry_dist)
    {
        /* there is an opponent near enough to worry about */
        DebugKick(printf("In RotToAvoidOpponent, avoiding opponent\n"));
        AngleDeg ball_to_targ = GetNormalizeAngleDeg((abs_dir - Mem->MyBodyAng()) - ball_ang);
        AngleDeg opp_to_targ = GetNormalizeAngleDeg((abs_dir - Mem->MyBodyAng()) - opp_ang);

        if (ball_to_targ * opp_to_targ < 0)
        { /* they're on opposite sides of the target */
            Mem->LogAction2(120, "RotToAvoidOpponent: CLOSEST");
            rot = TURN_CLOSEST;
        }
        else if (ball_to_targ < opp_to_targ)
        {
            Mem->LogAction2(120, "RotToAvoidOpponent: CW");
            rot = TURN_CW;
        }
        else
        {
            Mem->LogAction2(120, "RotToAvoidOpponent: CCW");
            rot = TURN_CCW;
        }
    }
    else
    {
        Mem->LogAction2(120, "RotToAvoidOpponent: no opponents close enough");
        rot = TURN_CLOSEST;
    }

    return rot;
}

/* Picks the rotation direction (CW or CCW) such that the ball has to travel
   the least distance */
TurnDir RotClosest(float abs_dir)
{
    AngleDeg cw_ang = (abs_dir - Mem->MyBodyAng()) - Mem->BallAngleFromBody();
    AngleDeg ccw_ang = Mem->BallAngleFromBody() - (abs_dir - Mem->MyBodyAng());
    DebugKick(printf("Test1: cw_ang: %f\tccw_ang: %f\n", cw_ang, ccw_ang));
    if (cw_ang < 0)
        cw_ang += 360;
    if (ccw_ang < 0)
        ccw_ang += 360;
    DebugKick(printf("Test2: cw_ang: %f\tccw_ang: %f\n", cw_ang, ccw_ang));
    if (cw_ang < ccw_ang)
        return TURN_CW;
    else
        return TURN_CCW;
}

/* Kick with direction ddir and power such that the ball moves distance ddist
   If ddist is too big, kick it as hard as possible.
   corrects for ball velocity
   Returns a command object partially filled out
   The distance is multiplied by the distFactor */
/* returns KickCommand.time = -1 for error */
TurnKickCommand dokick(AngleDeg ddir, float ddist, float distFactor,
                       Bool *pCanKickToTraj)
{
    float v0;
    float power;
    float kick_dir = ddir;
    TurnKickCommand com;
    com.time = -1;
    com.type = CMD_kick;
    com.turn_neck = FALSE;
    if (pCanKickToTraj)
        *pCanKickToTraj = TRUE;

    Mem->LogAction6(110, "dokick: %.2f (%.2f) at %.2f, rate %.5f",
                    ddist, distFactor, ddir, Mem->BallKickRate());

    v0 = ddist * distFactor;

    NormalizeAngleDeg(&ddir);

    DebugKick(printf(" dokick: ddir: %f\tv0: %f\n", ddir, v0));
    DebugKick(printf(" kickrate: %f\tmyang: %f\n", Mem->BallKickRate(), Mem->MyAng()));

    if (Mem->BallKickRate() == 0.0)
        my_error("dokick: Huh? BallKickRate is 0!");

    if (!Mem->BallVelocityValid())
    {
        DebugKick(printf("In dokick with velocity not valid. Assuming it's 0."));
        Mem->LogAction2(130, "dokick: assuming vel 0");
    }

    if (Mem->BallVelocityValid() &&
        Mem->BallAbsoluteVelocity() != 0)
    { /* correct for ball velocity */
        Vector tmp = Polar2Vector(v0, ddir) - Mem->BallRelativeToBodyVelocity();
        kick_dir = tmp.dir();
        power = tmp.mod() / Mem->BallKickRate();
        DebugKick(printf(" Correcting for ball velocity# vel.x: %f\tvel.y:%f\n",
                         Mem->BallRelativeVelocity().x, Mem->BallRelativeVelocity().y));
        DebugKick(printf(" New stuff# power: %f\t ddir:%f\n", power, kick_dir));
        Mem->LogAction4(130, "dokick: vel corr: %f at %f", power, kick_dir);
    }
    else
    {
        Mem->LogAction2(130, "dokick: vel is 0");
        power = v0 / Mem->BallKickRate();
    }

    if (power > Mem->SP_max_power)
    {
        DebugKick(printf("Trying to kick over SP_max_power! Correcting...\n"));
        // Mem->LogAction2(130, "dokick: trying to kick too hard, correcting");
        if (!Mem->BallVelocityValid() || Mem->BallAbsoluteVelocity() == 0)
        {
            power = Mem->SP_max_power;
            Mem->LogAction4(130, "dokick: max_pow corr: (stopped ball) %f at %f",
                            power, kick_dir);
        }
        else if (ddist == 0)
        {
            /* this is probably a stop_ball kick, but ddir is meaningless for
           the RayCircleIntersection below.
           Therefore, we just kick against the ball's velocity */
            kick_dir = GetNormalizeAngleDeg(Mem->BallRelativeToBodyHeading() + 180);
            power = Mem->SP_max_power;
            Mem->LogAction4(130, "dokick: max_pow corr: (ddist = 0) %f at %f",
                            power, kick_dir);
        }
        else
        {
            /* this is a ray circle intersection problem
               We want to kick in the right direction, but not as hard as desired */
            Vector sol1, sol2;
            int numSol;
            Vector vNewTraj;
            numSol =
                RayCircleIntersect(Ray(Mem->BallRelativeToBodyPosition(), ddir),
                                   Mem->SP_max_power * Mem->BallKickRate(),
                                   Mem->Global2RelativeToMyBody(Mem->BallPredictedPosition()),
                                   &sol1, &sol2);
            /* we want the solution that's furthest along the ray - that's sol2
           if there are two solution */
            if (numSol == 0)
            {
                /* we can't kick ball to desired trajectory, so let's just slow it
                   down by kicking directly against velocity
                   It might be better to return an error and have the caller decide
                   what to do, but this works pretty well */
                DebugKick(printf("Can't kick ball to right trajectory!\n"));
                power = Mem->SP_max_power;
                kick_dir = Mem->BallRelativeToBodyHeading() + 180;
                if (pCanKickToTraj)
                    *pCanKickToTraj = FALSE;
                Mem->LogAction4(130, "dokick: max_pow corr: can't get to trajectory! %.2f at %.2f",
                                power, kick_dir);
            }
            else
            {
                if (numSol == 1)
                {
                    vNewTraj = sol1;
                }
                else if (numSol == 2)
                {
                    /* we want the solution closer to the target point we wanted */
                    Vector targ = Mem->BallRelativeToBodyPosition() + Polar2Vector(v0, ddir);
                    Mem->LogAction8(130, "dokick: 2 solutions (%.2f, %.2f) (%.2f, %.2f) targ (%.2f, %.2f)",
                                    sol1.x, sol1.y, sol2.x, sol2.y, targ.x, targ.y);
                    if (sol1.dist2(targ) < sol2.dist(targ))
                    {
                        Mem->LogAction4(140, "Picked sol1: %.2f < %.2f",
                                        sol1.dist2(targ), sol2.dist2(targ));
                        vNewTraj = sol1;
                    }
                    else
                    {
                        Mem->LogAction4(140, "Picked sol2: %.2f > %.2f",
                                        sol1.dist2(targ), sol2.dist2(targ));
                        vNewTraj = sol2;
                    }
                }
                else
                {
                    my_error("dokick: How many solutions to RayCircleIntersection? %d", numSol);
                }

                Vector vNewKick = vNewTraj -
                                  Mem->Global2RelativeToMyBody(Mem->BallPredictedPosition());
                power = vNewKick.mod() / Mem->BallKickRate();
                kick_dir = vNewKick.dir();

                DebugKick(printf(" Correcting for ball velocity AND max power# vel.x: %f\tvel.y:%f\n",
                                 Mem->BallRelativeVelocity().x, Mem->BallRelativeVelocity().y));
                DebugKick(printf("New stuff# power: %f\t dir:%f\n", power, kick_dir));
                Mem->LogAction4(130, "dokick: max_pow corr: %f at %f", power, kick_dir);
            }
        }
    }

    power = Min(Round(power, -2), Mem->SP_max_power);
    kick_dir = Round(kick_dir, -2);
    NormalizeAngleDeg(&kick_dir);
    DebugKick(printf("kicking with power: %f at dir: %f\n", power, kick_dir));
    com.time = Mem->CurrentTime;
    com.power = power;
    com.angle = kick_dir;
    return com;
}

/* all the turnball reasoning is done in relative (to body) coordinates
   call this before calling dokick to correct for movement of the player */
void PlayerMovementCorrection(AngleDeg *pdir, float *pdist)
{
    Mem->LogAction4(110, "Player movement correction: before dist: %f at dir: %f",
                    *pdist, *pdir);
    DebugKick(printf("Before corr- dir: %f\t dist: %f\n", *pdir, *pdist));
    Vector vRelBallTarg = Mem->BallRelativeToBodyPosition() +
                          Polar2Vector(*pdist, *pdir);
    Vector vNewRelTarg = Mem->MyPredictedPosition() +
                         vRelBallTarg - Mem->MyPos();
    Vector vNewRelTraj = vNewRelTarg - Mem->BallRelativeToBodyPosition();
    *pdir = vNewRelTraj.dir();
    *pdist = vNewRelTraj.mod();
    DebugKick(printf("After corr- dir: %f\t dist: %f\n", *pdir, *pdist));
    Mem->LogAction4(110, "Player movement correction: after dist: %f at dir: %f",
                    *pdist, *pdir);
}

/* kick ball so that it is at angle ddir and dist EndDist
   If you have to kick around the player, kick rotway(clockwise or counter-)
   */
/* we always reason about the right trajectory for the ball leave velocity
   correction for dokick */
KickToRes turnball_kick(AngleDeg target_dir, TurnDir rotate,
                        Bool StopBall, TurnKickCommand *pCom,
                        float EndDist, float closeMarg, float kickFac)
{
    float dir;
    float dist;
    Vector btraj;

    pCom->time = -1;
    pCom->turn_neck = FALSE;

    DebugKick(printf("\nat turnball_kick: target_dir: %f\n", target_dir));
    Mem->LogAction4(60, "Turnball_kick: targ_dir: %.1f  dir: %d", target_dir, (int)rotate);

    NormalizeAngleDeg(&target_dir);

    // DebugKick(printf("HERE Time: %d\n", Mem->CurrentTime.t));
    /* the pos valid is not really right - if we are turning the ball and didn't
       actually see it last time, then there's a problem */
    if (!Mem->BallPositionValid() || !Mem->BallKickable())
    {
        Mem->LogAction2(90, "turnball_kick: lost the ball");
        return KT_LostBall;
    }

    /* if the velocity isn's valid, turn to face ball */
    if (!Mem->BallVelocityValid())
    {
        float ball_ang_from_body = Mem->BallAngleFromBody(); /* for efficiency */
        Mem->LogAction2(90, "turnball_kick: vel is not valid, looking at it");
        DebugKick(printf("turning to face ball\n"));
        if (Mem->CanSeeBallWithNeck())
        {
            pCom->time = Mem->CurrentTime;
            pCom->type = CMD_kick;
            pCom->angle = ball_ang_from_body + 180;
            pCom->power = Mem->CP_stop_ball_power;

            pCom->turn_neck = TRUE;
            pCom->turn_neck_angle = Mem->LimitTurnNeckAngle(Mem->BallAngleFromNeck());
        }
        else
        {
            /* turn body to face ball, and turn neck to straight ahead */
            pCom->time = Mem->CurrentTime;
            pCom->type = CMD_turn;
            pCom->turn_neck = TRUE;
            if (fabs(ball_ang_from_body) > Mem->MaxEffectiveTurn())
            {
                /* out body can't get to where we want to go */
                pCom->angle = 180; /* get our maximum effective turn */
                pCom->turn_neck_angle = ball_ang_from_body -
                                        signf(ball_ang_from_body) * Mem->MaxEffectiveTurn();
            }
            else
            {
                pCom->angle = ball_ang_from_body;
                pCom->turn_neck_angle = -Mem->MyNeckRelAng();
            }
        }

        return KT_TurnedToBall;
    }

    DebugKick(printf(" ball.dist: %f\t.dir: %f\n",
                     Mem->BallDistance(), Mem->BallAngle()));
    DebugKick(printf(" HERE ball.vel.x: %f\t.y: %f\tmod: %f\n",
                     Mem->BallRelativeVelocity().x, Mem->BallRelativeVelocity().y,
                     Mem->BallSpeed()));
    DebugKick(printf(" ball.rpos.x: %f\t.y: %f\n",
                     Mem->BallRelativePosition().x, Mem->BallRelativePosition().y));
    DebugKick(printf(" target_dir: %f\n", target_dir));

    if (fabs(GetNormalizeAngleDeg(target_dir - Mem->BallAngleFromBody())) < Mem->CP_KickTo_err)
    {
        /* Do something to indicate we are done */
        if (!StopBall || Mem->BallSpeed() < Mem->CP_max_ignore_vel)
            return KT_DidNothing;
        Mem->LogAction2(90, "turnball_kick: we're there, stopping the ball");
        DebugKick(printf("  Stop ball kick\n"));
        dir = 0;
        dist = 0;
        PlayerMovementCorrection(&dir, &dist);
        *pCom = dokick(dir, dist, 1.0);
        pCom->turn_neck = FALSE;
        return KT_Success;
    }

    if (rotate == TURN_AVOID)
    {
        rotate = RotToAvoidOpponent(target_dir + Mem->MyBodyAng());
    }

    if (rotate == TURN_CLOSEST)
    {
        rotate = RotClosest(target_dir + Mem->MyBodyAng());
    }

    if (is_straight_kick(target_dir, EndDist, closeMarg))
    {
        float pow;

        btraj = Polar2Vector(EndDist, target_dir) - Mem->BallRelativeToBodyPosition();
        dir = btraj.dir();
        dist = btraj.mod();

        /* now we're goign to do some distance twiddling to get the ball to
           get to the right angle and stop */
        pow = dist / Mem->BallKickRate();
        pow = Min(pow, Mem->CP_max_turn_kick_pow);
        dist = pow * Mem->BallKickRate();

        Mem->LogAction4(90, "turnball_kick: striaght kick: dist %f at %f", dist, dir);
        DebugKick(printf("  Straight kick# dir: %f dist: %f\n", dir, dist));
        PlayerMovementCorrection(&dir, &dist);
        *pCom = dokick(dir, dist, 1.0);
        pCom->turn_neck = FALSE;
    }
    else if (Mem->BallDistance() < closeMarg)
    {

        /* ball is too close to do a tangent kick, so do a kick at 90 degrees */
        dir = ((int)rotate) * (-90) + Mem->BallAngleFromBody();
        dist = 2.0 * sqrt(Sqr(Mem->CP_opt_ctrl_dist) - Sqr(Mem->BallDistance()));
        Mem->LogAction2(90, "turnball_kick: 90 deg kick");
        DebugKick(printf("  Close kick# dir: %f dist: %f\n", dir, dist));
        PlayerMovementCorrection(&dir, &dist);
        *pCom = dokick(dir, dist, kickFac);
        pCom->turn_neck = FALSE;
    }
    else
    {

        /* do a turning kick */
        /* we make a circle around the player of radius closeMarg
           and calculate the trajectory that goes in the right direction and is
           tangent to the circle */
        dir = 180 + Mem->BallAngleFromBody() + ((int)rotate) * ASin(closeMarg / Mem->BallDistance());
        DebugKick(printf(" ball dist: %f\tclosest_margin: %f\n",
                         Mem->BallDistance(), closeMarg));
        dist = sqrt(Sqr(Mem->BallDistance()) - Sqr(closeMarg));
        dist +=
            sqrt(Sqr(Mem->CP_opt_ctrl_dist) - Sqr(closeMarg));
        DebugKick(printf("  Turning ball# dir: %f dist: %f\n", dir, dist));
        Mem->LogAction2(90, "turnball_kick: turning kick");
        PlayerMovementCorrection(&dir, &dist);
        *pCom = dokick(dir, dist, kickFac);
        pCom->turn_neck = FALSE;
    }

    return KT_DidKick;
}

KickToRes TurnballTo(AngleDeg rel_dir, TurnDir rotate)
{
    TurnKickCommand com;
    KickToRes res = turnball_kick(rel_dir, rotate, TRUE, &com);
    if (res == KT_Success || res == KT_DidKick || res == KT_TurnedToBall)
        DoTurnKickCommand(com);

    return res;
}

/*******************************************************************************************/

/* if we kick the ball as hard as possible in the right direction, will it
   be a collision with the player? */
/* all we have to do is look at the predicted ball position with that kick and
   see if it is within the player's radius of the player */
int is_hard_kick_coll(float abs_dir, TurnKickCommand *pcom,
                      Vector *pvPredBall, float targ_vel, Bool *pCanKickToTraj = NULL)
{
    DebugKick2(cout << "Is hard_kick_coll here" << endl);
    *pcom = dokick(abs_dir - Mem->MyBodyAng(), targ_vel, 1.0, pCanKickToTraj);
    *pvPredBall = Mem->BallPredictedPosition(1, pcom->power, pcom->angle);
    DebugKick(cout << "IsColl: PredBall: " << *pvPredBall
                   << "\tMyPos: " << Mem->MyPredictedPosition() << endl);
    DebugKick(cout << "diff: " << (*pvPredBall - Mem->MyPredictedPosition()).mod()
                   << "\tmarg: " << Mem->CP_hard_kick_margin << endl);
    return (*pvPredBall - Mem->MyPredictedPosition()).mod() <=
           Mem->SP_player_size + Mem->CP_hard_kick_dist_buffer;
}

/* Used when we decide a kick in the right direction woudl be a collision,
   so we need to turnball to kick the ball more to the side of us,
   OR, when in KM_HardestKick, moving the ball back for a kick */
TurnKickCommand hard_kick_turnball(float abs_dir, TurnDir rot, Bool StopBall = FALSE)
{
    TurnKickCommand com;
    // TurnDir rot = RotToAvoidOpponent(abs_dir);
    /* SMURF - should this have a larger dokick_factor? */
    Mem->LogAction3(70, "hard_kick_turnball: %f", abs_dir);
    DebugKick2(cout << "Doing a hard_kick_turnball" << endl);
    KickToRes res = turnball_kick(abs_dir - Mem->MyBodyAng(), rot, StopBall,
                                  &com, Mem->CP_hardest_kick_ball_dist,
                                  Mem->CP_closest_margin, Mem->CP_dokick_factor);
    if (res == KT_DidNothing || res == KT_LostBall)
        my_error("	hard_kick_turnball: Something weird happened: %d", res);
    return com;
}

TurnKickCommand kick_hard_moderate(AngleDeg abs_dir, float targ_vel, TurnDir rot)
{
    /* Here's the idea:
         See if one strong kick will be a collision (if it is, turnball)
         o.w. manipulate our kick so that we get another one next cycle
         (but make sure the final velocity is higher)
         it that makes it a collision, just do the strong kick */
    /* there is no reasning about turning the ball backward to get max power.
       See KM_HardestKick for that */
    Mem->LogAction2(60, "kick_hard_moderate");
    DebugKick2(cout << endl
                    << "Time: " << Mem->CurrentTime.t
                    << "\tkick_hard_moderate called" << endl);
    TurnKickCommand kickCom;
    TurnKickCommand HKCommand;
    Vector vPredBall;
    Bool CanKickToTraj;
    if (is_hard_kick_coll(abs_dir, &kickCom, &vPredBall, targ_vel, &CanKickToTraj))
    {
        DebugKick2(cout << " Moderate: collision, still thinking about it" << endl);
        if (Mem->BallDistance() < Mem->SP_player_size + Mem->CP_hard_kick_dist_buffer)
        {
            Mem->LogAction2(70, "kick_hard_moderate: ball is too close to plan for collision ");
            return hard_kick_turnball(abs_dir, rot);
        }
        else
        {
            if (CanKickToTraj)
            {
                Mem->LogAction2(70, "kick_hard_moderate: planning for collision ");

                /* figure out what kick it would be to just go to the edge of the player */
                Vector sol1, sol2;
                int numSol =
                    RayCircleIntersect(Ray(Mem->BallAbsolutePosition(),
                                           vPredBall - Mem->BallAbsolutePosition()),
                                       Mem->SP_player_size + Mem->CP_hard_kick_dist_buffer,
                                       Mem->MyPredictedPosition(),
                                       &sol1, &sol2);
                if (numSol == 0)
                {
                    my_error("Woah! No solution to kicking a ball less hard!");
                    return hard_kick_turnball(abs_dir, rot);
                }
                Vector vBallTarg = sol1;
                /* SMURF: is there any case where we wouldn't want to do this kick? */
                DebugKick2(cout << "Trying to kick to edge of player" << endl);
                Vector vBallPath = vBallTarg.Global2Relative(Mem->BallAbsolutePosition(), Mem->MyBodyAng());
                HKCommand = dokick(vBallPath.dir(), vBallPath.mod());
                return HKCommand;
            }
            else
            {
                /* can't kick to trajectory and it's a collision.
                   That's fine! Just do the kick */
                Mem->LogAction2(70, "kick_hard_moderate: collision, but that's okay!");
                return kickCom;
            }
        }
    }
    else
    {
        DebugKick(cout << " MyPos: " << Mem->MyPos() << endl);
        DebugKick(cout << " vPredBall: " << vPredBall << endl);
        Vector vPredBallRel =
            vPredBall.Global2Relative(Mem->MyPredictedPosition(), Mem->MyBodyAng());
        DebugKick(cout << " vPredBallRel: " << vPredBallRel << endl);
        if (vPredBallRel.mod() < Mem->SP_kickable_area - Mem->CP_hard_kick_dist_buffer ||
            !Mem->BallWillBeKickable())
        {
            /* we're going to get another kick next time or this will be our last
             kick anyway - do the strong one! */
            Mem->LogAction2(70, "kick_hard_moderate: strong with another or last kick");
            DebugKick2(cout << " Moderate: strong with another kick or last kick!" << endl);
            return kickCom;
        }
        DebugKick2(cout << " Moderate: deciding whether to do one or two kicks" << endl);
        /* we're goign to set vBall to be the relative vector (to new pos)
         to the position that will give us another kick next time */
        float oneKickVel = Min(Mem->SP_ball_speed_max,
                               (vPredBall - Mem->BallAbsolutePosition()).mod());
        float twoKickVel = 0.0;
        Vector sol1, sol2;
        int numSol;
        Vector vBallTarg;
        numSol =
            RayCircleIntersect(Ray(Mem->BallAbsolutePosition(),
                                   vPredBall - Mem->BallAbsolutePosition()),
                               Mem->SP_kickable_area - Mem->CP_hard_kick_dist_buffer,
                               Mem->MyPredictedPosition(),
                               &sol1, &sol2);
        /* we want the solution that's furthest along the ray - that's sol2
         if there are two solution */
        if (numSol != 0)
        {

            if (numSol == 2)
                vBallTarg = sol2;
            else
                vBallTarg = sol1;

            /* see if the first of the two kicks is a coll */
            if ((vBallTarg - Mem->MyPredictedPosition()).mod() >=
                Mem->SP_player_size + Mem->CP_hard_kick_dist_buffer)
            {
                /* we can do it without collision */
                /* now see if this is actually goign to be better */
                vBallTarg =
                    vBallTarg.Global2Relative(Mem->BallAbsolutePosition(), Mem->MyBodyAng());
                float kickrate =
                    Mem->GetBall()->calc_kick_rate(vBallTarg.mod(), vBallTarg.dir());
                /* the first kick */
                // DebugKick( cout << "vBallTarg: " << vBallTarg << "\tBallAbsPos: " <<
                //	 Mem->BallAbsolutePosition() << endl );
                twoKickVel = (vBallTarg - Mem->BallRelativeToBodyPosition()).mod() *
                             Mem->SP_ball_decay;
                DebugKick(printf("  twoKickVel: first kick: %f\n", twoKickVel));
                /* the second kick */
                twoKickVel = Min(Mem->SP_ball_speed_max,
                                 twoKickVel + kickrate * Mem->SP_max_power);
                DebugKick(printf("  oneKickVel: %f\ttwoKickVel: %f\n", oneKickVel, twoKickVel));
            }
            else
                my_error("kick_hard_moderate- no ray intersection?");
        }

        /* remember if twoKick is a collision, then it's velocity will be 0 */
        if (numSol == 0 || oneKickVel >= twoKickVel)
        {
            /* do the one hard kick */
            Mem->LogAction2(70, "kick_hard_moderate: doing one hard kick");
            DebugKick2(cout << " Moderate- Doing one hard kick" << endl);
            return kickCom;
        }
        else
        {
            /* do the weaker kick */
            Mem->LogAction2(70, "kick_hard_moderate: doing first of two kicks");
            DebugKick2(cout << " Moderate- doing first of two kicks" << endl);
            // DebugKick(printf(" Predicted distance: %f\n",
            //		 (vBallTarg - Mem->MyPredictedPosition()).mod() ));
            // DebugKick(cout << " BallCurrPos: " << Mem->BallAbsolutePosition() << endl);
            HKCommand = dokick(vBallTarg.dir(), vBallTarg.mod());
            DebugKick(cout << " KickTraj: " << vBallTarg << endl);
            DebugKick(cout << " PredPos: " << Mem->MyPredictedPosition() << endl);
            return HKCommand;
        }
    }
}

TurnDir KickRotationDirectionAbs(AngleDeg abs_ang, TurnDir rot)
{
    if (rot == TURN_AVOID)
        rot = RotToAvoidOpponent(abs_ang);
    if (rot == TURN_CLOSEST)
        rot = RotClosest(abs_ang);
    if (rot != TURN_CW && rot != TURN_CCW)
        my_error("KickRotationDirection: bad TurnDir");

    return rot;
}

/* see above for description of KM_Moderate */
/* KM_HardestKick: Moves the ball to the side of us (relative to the kick
   direction) then uses KM_Moderate */
/* KM_Quickly, KM_QuickestRelease: get rid of the ball as fast as possible,
   will turnball if neccesary. KM_Quickly will turn to face the ball if it can
   so that we get a harder kick */
/* returns 1 if a kick actually done */
int smart_kick_hard_abs(float abs_dir, KickMode mode, float targ_vel,
                        TurnDir rot)
{
    TurnKickCommand HKCommand;
    HKCommand.time = -1;
    HKCommand.turn_neck = FALSE;

    DebugKick(printf("\nsmart_kick_hard: Time: %d\n", Mem->CurrentTime.t));

    Mem->LogAction4(50, "smart_kick_hard_abs: angle = %.1f, mode = %d", abs_dir, mode);

    if (!Mem->BallPositionValid())
    {
        my_error("smart_kick_hard called with ball position not valid");
        return 0;
    }

    if (!Mem->BallKickable())
    {
        my_error("smart_kick_hard called with ball not kickable!");
        return 0;
    }

    /* if the velocity isn's valid, turn to face ball */

    if (!Mem->BallVelocityValid())
    {
        float ball_ang_from_body = Mem->BallAngleFromBody(); /* for efficiency */
        Mem->LogAction2(60, "smart_kick_hard_abs: turning to face ball");
        DebugKick2(printf("smart_kcik_hard: turning to face ball\n"));
        if (Mem->CanSeeBallWithNeck())
        {
            Mem->LogAction2(70, "smart_kick_hard_abs: just turniong neck");
            HKCommand.time = Mem->CurrentTime;
            HKCommand.type = CMD_kick;
            HKCommand.angle = ball_ang_from_body + 180;
            HKCommand.power = Mem->CP_stop_ball_power;

            HKCommand.turn_neck = TRUE;
            HKCommand.turn_neck_angle = Mem->LimitTurnNeckAngle(Mem->BallAngleFromNeck());
        }
        else
        {
            /* turn body to face ball, and turn neck to straight ahead */
            Mem->LogAction2(70, "smart_kick_hard_abs: turning neck and body");
            HKCommand.time = Mem->CurrentTime;
            HKCommand.type = CMD_turn;
            HKCommand.turn_neck = TRUE;
            if (fabs(ball_ang_from_body) > Mem->MaxEffectiveTurn())
            {
                /* out body can't get to where we want to go */
                HKCommand.angle = 180; /* get our maximum effective turn */
                HKCommand.turn_neck_angle = ball_ang_from_body -
                                            signf(ball_ang_from_body) * Mem->MaxEffectiveTurn();
            }
            else
            {
                HKCommand.angle = ball_ang_from_body;
                HKCommand.turn_neck_angle = -Mem->MyNeckRelAng();
            }
        }
        return DoTurnKickCommand(HKCommand);
    }

    rot = KickRotationDirectionAbs(abs_dir, rot);

#ifdef NEVER
    /* With the new kick code, the ball goes through the player,
       so we don't usually need to do a turn-ball as part of a shot */
    if (mode <= KM_Moderate &&
        Mem->IsPointBehind(Mem->BallAbsolutePosition(), abs_dir))
    {
        /* see if we need to rotate one way */
        DebugKick(cout << "smart_kick_hard: decign if rotating " << rot << endl);
        /* now decide if ball is on the wrong side	 */
        if (Mem->IsPointBehind(Mem->BallAbsolutePosition(), abs_dir + ((int)rot) * 90))
        {
            /* start rotating right way */
            DebugKick(cout << "smart_kick_hard: special turnball to avoid opp" << endl);
            TurnKickCommand com;
            KickToRes res = turnball_kick(abs_dir - Mem->MyAng(),
                                          // abs_dir + ((int)rot)*90 - Mem->MyAng(),
                                          rot, FALSE, &com);
            if (res == KT_DidNothing || res == KT_LostBall)
                my_error("smart_kick_hard: special turnball; turnkick failed, res: %d", res);
            else
                return DoTurnKickCommand(com);
        }
    }
#endif

    switch (mode)
    {
    case KM_None:
        my_error("KM_None is not a valid kick mode for smart_kick_hard!");
        break;

    case KM_HardestKick:
    {
        if (Mem->CurrentTime - 1 != Mem->HKTime)
        {
            DebugKick(printf("First in a chain of calls\n"));
            /* this is the first in a chain of calls */
            Mem->HKStep = 0;
            /* decide which way to rotate ball */
            Mem->HKrot = rot;
            /*AngleDeg ang = Mem->BallAngle() + Mem->MyAng() - abs_dir;
            NormalizeAngleDeg(&ang);
            if (ang >= 0)
          Mem->HKrot = TURN_CW;
            else
          Mem->HKrot = TURN_CCW;*/

            /* see if we need to turn */
            AngleDeg target_dir = abs_dir - Mem->MyBodyAng() + ((int)Mem->HKrot) * 90;
            if (fabs(target_dir) > Mem->CP_max_hard_kick_angle_err)
            {
                Mem->LogAction2(70, "smart_kick_hard_abs: hardest turning for power");
                DebugKick(printf("turning to get hard kick\n"));
                Mem->HKTime = Mem->CurrentTime;
                Mem->HKStepNext = 0;
                HKCommand.type = CMD_turn;
                HKCommand.angle = target_dir;
                HKCommand.time = Mem->CurrentTime;
                HKCommand.turn_neck = FALSE;
                break; // used to be: return DoTurnKickCommand(HKCommand);
            }
        }
        else
            Mem->HKStep = Mem->HKStepNext;

        /* AngleDeg turn_target = abs_dir + 180 +
          ((int)Mem->HKrot)*Mem->CP_hardest_kick_angle_disp; */
        AngleDeg turn_target = abs_dir +
                               ((int)Mem->HKrot) * (Mem->CP_hardest_kick_ball_ang);
        NormalizeAngleDeg(&turn_target);
        if (Mem->HKStep == 0 &&
            fabs(GetNormalizeAngleDeg(turn_target - Mem->MyBodyAng() - Mem->BallAngleFromBody())) >
                Mem->CP_KickTo_err)
        {
            /* on step 0, we turn ball to back of us */
            Mem->LogAction2(70, "smart_kick_hard_abs: hardest turning ball backwards");
            HKCommand = hard_kick_turnball(turn_target, TURN_CLOSEST, TRUE); // TRUE to stop ball
            Mem->HKStepNext = 0;
        }
        if (Mem->HKStep == 1 ||
            (Mem->HKStep == 0 && HKCommand.time != Mem->CurrentTime))
        {
            /* on step 1, we call the moderate code */
            /* or if step 0 had turnball problems, we'll go ahead and drop
           to this code */
            Mem->LogAction2(70, "smart_kick_hard_abs: calling kick_hard_moderate");
            HKCommand = kick_hard_moderate(abs_dir, targ_vel, rot);
            Mem->HKStepNext = 1;
        }
        if (Mem->HKStep != 0 && Mem->HKStep != 1)
            my_error("HKstep is not in a valid state: %d", Mem->HKStep);

        Mem->HKTime = Mem->CurrentTime;
        break; // used to be: return DoTurnKickCommand(HKCommand);
    }

    case KM_Hard:
        /* SMURF */
        my_error("KM_Hard not implemented yet!");
        // dump_core("blah");
        break;

    case KM_Moderate:
        HKCommand = kick_hard_moderate(abs_dir, targ_vel, rot);
        break;

    case KM_Quickly:
    case KM_QuickestRelease:
    {
        Mem->HKStep = Mem->HKStepNext = -1;
        /* see if the hardest kick in the direction we want will be collision */
        /* use our dokick function to correct for vel */
        TurnKickCommand kickCom;
        Vector vPredBall;
        Bool CanKickToTraj;
        if (!is_hard_kick_coll(abs_dir, &kickCom, &vPredBall, targ_vel, &CanKickToTraj))
        {
            /* there is no collsion! */
            if (mode == KM_Quickly &&
                Mem->BallWillBeKickable() &&
                fabs(Mem->BallAngleFromBody()) > Mem->CP_max_hard_kick_angle_err)
            {
                /* In KM_Quickly mode, we can turn to ball and get it next time */
                Mem->LogAction2(70, "smart_kick_hard_abs: km_quickly: turning for power");
                HKCommand.time = Mem->CurrentTime;
                HKCommand.type = CMD_turn;
                HKCommand.angle = Mem->AngleToFromBody(Mem->BallPredictedPosition());
            }
            else
                HKCommand = kickCom;
            break; // used to be:return DoTurnKickCommand(HKCommand);
        }
        else
        {
            /* do turnball */
            if (!CanKickToTraj)
            {
                /* can't kick to trajectory, so just do the collision kick */
                Mem->LogAction2(70, "smart_kick_hard_abs: km_quick: can't kick to traj, taking the collision");
                HKCommand = kickCom;
            }
            else
            {
                Mem->LogAction2(70, "smart_kick_hard_abs: km_quick: turnball");
                HKCommand = hard_kick_turnball(abs_dir, rot);
            }
            break; // used to be:return DoTurnKickCommand(hard_kick_turnball(abs_dir, rot));
        }
    }
        my_error("How did I get to end of KM_Quickly/KM_Quickest_Release");
        break;

    default:
        my_error("Invalid/Unimplemented kick mode passed to smart_kick_hard");
        break;
    }

    /* now we'll add a turn_neck into this to look at where the ball will be
       If a turn_neck is already there, we won't do this */
    if (!HKCommand.turn_neck)
    {
        if (HKCommand.time != Mem->CurrentTime)
            my_error("Adding a turn_neck to an invalid HKCommand");
        Mem->LogAction2(70, "smart_kick_hard_abs: doing a turn_neck");
        Vector pred_ball_pos;
        if (HKCommand.type == CMD_kick)
            pred_ball_pos = Mem->BallPredictedPosition(1, HKCommand.power, HKCommand.angle);
        else
            pred_ball_pos = Mem->BallPredictedPosition(1);
        float pred_ball_dir = (pred_ball_pos - Mem->MyPredictedPosition()).dir();
        pred_ball_dir -= Mem->MyNeckGlobalAng();
        pred_ball_dir -= (HKCommand.type == CMD_turn ? HKCommand.angle : 0);

        HKCommand.turn_neck = TRUE;
        HKCommand.turn_neck_angle = Mem->LimitTurnNeckAngle(GetNormalizeAngleDeg(pred_ball_dir));
    }

    return DoTurnKickCommand(HKCommand);
}

/* good passes require that the ball is not moving too quickly when it reaches
   the intended recipient, so this cover function helps acheive that */
int smart_pass(Vector pt, float targ_vel_at_pt, KickMode mode, TurnDir rot)
{
    Mem->LogAction2(60, "smart_pass: doing a turn_neck");
    return smart_kick_hard(Mem->AngleToFromBody(pt), mode,
                           Mem->VelAtPt2VelAtFoot(pt, targ_vel_at_pt),
                           rot);
}

/* should be used when we have free kicks and stuff */
/* kick_ang: absolute angle for the direction that the ball will be kicked */
/* returns whether we are in the right place or not */
Bool go_to_static_ball(float kick_ang)
{
    Line l;

    Mem->LogAction5(50, "go_to_static_ball: ball at (%.1f, %.1f) for kick angle %.2f",
                    Mem->BallX(), Mem->BallY(), kick_ang);

    if (!Mem->BallPositionValid())
        my_error("go_to_static_ball: lost ball");

    /* we can be pretty tolerant of angle errors, but distance errors start to matter
       real quickly */

    Vector targ_pt = Mem->BallAbsolutePosition() +
                     Polar2Vector(Mem->CP_hardest_kick_ball_dist, kick_ang);

    /* we want to try and face the ball at all times */
    turn_neck(Mem->LimitTurnNeckAngle(Mem->BallAngleFromNeck()));

    if (Mem->MySpeed() == 0 &&
        Mem->DistanceTo(targ_pt) <= Mem->CP_static_kick_dist_err)
    {
        return TRUE;
    }

    /* if we are real far from the ball, just use the regular go_to_point */
    if (Mem->BallDistance() > 2 * Mem->SP_kickable_area)
    {
        Mem->LogAction2(60, "go_to_static_ball: far away, using go_to_point");
        if (go_to_point(targ_pt, 0 /* no buffer */, Min(Mem->SP_stamina_inc, Mem->SP_max_power)) != AQ_ActionQueued)
            my_error("go_to_static_ball: far away, why didn't go_to_point do anything?");
        return FALSE;
    }

    /* make sure that we go around the ball */
    l.LineFromTwoPoints(Mem->MyPos(), targ_pt);
    Vector proj_ball_pt = l.ProjectPoint(Mem->BallAbsolutePosition());
    if (proj_ball_pt.dist(Mem->BallAbsolutePosition()) <=
            Mem->SP_player_size + Mem->SP_ball_size + Mem->CP_collision_buffer &&
        l.InBetween(proj_ball_pt, Mem->MyPos(), targ_pt))
    {
        /* we'll do a 90 degree dodge -we always go right */
        Vector dodge_pt = Mem->MyPos() +
                          Polar2Vector(Mem->SP_player_size, Mem->BallAngleFromBody() + Mem->MyBodyAng() + 90);
        Mem->LogAction2(60, "go_to_static_ball: dodging the ball");
        if (go_to_point(dodge_pt, 0 /* no buffer */, Min(Mem->SP_stamina_inc, Mem->SP_max_power)) != AQ_ActionQueued)
            my_error("go_to_static_ball: dodging, why didn't go_to_point do anything?");
        return FALSE;
    }

    /* now we need to get to the target_point */
    /* first see if we need to turn */
    l.LineFromRay(Mem->MyPos(), Mem->MyBodyAng());
    float ang = Mem->AngleToFromBody(targ_pt);
    if (fabs(ang) > 90 ||
        (l.dist(targ_pt) > Mem->CP_static_kick_dist_err &&
         ang > Mem->CP_static_kick_ang_err))
    {
        Mem->LogAction2(60, "go_to_static_ball: turning to target_point");
        turn(Mem->AngleToFromBody(targ_pt));
        return FALSE;
    }

    /* now calculate the speed we should be going to land right on the point */
    float targ_speed =
        SolveForFirstTermInfGeomSeries(Mem->SP_player_decay, Mem->DistanceTo(targ_pt));
    float dash_pow =
        MinMax(-Mem->SP_stamina_inc / 2,
               (targ_speed - Mem->MySpeed()) / Mem->SP_dash_power_rate,
               Mem->SP_stamina_inc);
    Mem->LogAction5(60, "go_to_static_ball: targ_speed: %.2f\tMySpeed: %.2f\tdash_pow: %.2f",
                    targ_speed, Mem->MySpeed(), dash_pow);
    if (fabs(dash_pow) > 1)
    {
        dash(dash_pow);
        return FALSE;
    }

    return TRUE;
}

/* -*- Mode: C++ -*- */

/* dribble.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

#ifdef DEBUG_OUTPUT
#define DebugDrib(x)
#define DebugDrib2(x)
#else
#define DebugDrib(x)
#define DebugDrib2(x)
#endif

void PatsTest_dribble()
{
    static Vector pos(10, 10);
    static AngleDeg drib_ang = 60;

    if (!Mem->MyConf())
    {
        Mem->LogAction2(10, "Don't know where I am, scanning field");
        scan_field_with_body();
        return;
    }

    if (!Mem->BallPositionValid())
    {
        Mem->LogAction2(10, "Lost ball, scanning field");
        return;
    }

    if (Mem->CurrentTime.t % 40 == 1)
        drib_ang = -drib_ang;

    if (!Mem->BallVelocityValid())
    {
        if (Mem->BallKickable())
            stop_ball();
        else
            get_ball();
    }
    else
    {
        /*printf("\nTime: %d\n", Mem->CurrentTime.t);*/
        DribbleRes res = DribbleTo(pos, Mem->CP_dribble_dash_pow,
                                   1, drib_ang, DM_Lazy);
        if (res == DR_GotThere)
        {
            Mem->LogAction2(10, "got to dribble target, reversing");
            stop_ball();
            pos = -pos;
        }
        else if (res == DR_LostBall)
        {
            get_ball();
        }
    }

    change_view(VW_Narrow);
}

/* Here's the basic strategy for dribbling:
   Alternate kicks and dashs so that the ball stays close to us.
   Using predicted ball positions and velocity corrected kicks, then given
   an angle we want to dribble at, we can choose a dash or kick that is most
   appropraite */

/* Checks to see if a dash would keep us in control of ball or if it would
   help us catch up to ball. If not, it calls dribble_kick */
TurnKickCommand dribble_dash(Vector vEndPos, float max_pow,
                             AngleDeg drib_ang, DribbleMode mode)
{
    Mem->LogAction2(70, "dribble: trying to do a dash");
    DebugDrib(cout << "in dribble_dash" << endl);
    TurnKickCommand com;
    com.time = -1;
    com.turn_neck = FALSE;

    if (!Mem->BallWillBeKickable(1, max_pow, Mem->CP_kickable_buffer))
    {
        if (Mem->BallKickable())
        {
            DebugDrib(printf("Dribble: Ball is not goign to be kickable, so we'll kick it.\n"));
            Mem->LogAction2(80, "dribble: tried dash, but will lose ball, so kicking");
            return dribble_kick(vEndPos, max_pow, drib_ang, mode);
        }
        if (!Mem->BallWillBeKickable(2, max_pow, Mem->CP_kickable_buffer) &&
            !Mem->BallWillBeKickable(3, max_pow, Mem->CP_kickable_buffer))
        {
            DebugDrib(printf("Dribble: looks bad for dash, probably need to go to intercept"));
            Mem->LogAction2(70, "dribble: lost ball it seems");
            return com;
        }
        else
            DebugDrib(printf("Dribble: It's okay, ball will be kickable soon\n"));
    }

    if (Mem->WillDashBeCollision(max_pow, Mem->CP_collision_buffer))
    {
        DebugDrib(printf("Dribble: dash would be a collision, so we'll kick instead\n"));
        if (Mem->BallKickable())
        {
            Mem->LogAction2(70, "dribble: dash would collide, so we kick");
            return dribble_kick(vEndPos, max_pow, drib_ang, mode);
        }
        else
        {
            /* we'll fall through to a dash-
           that'll at least get us close to the ball */
            DebugDrib(printf("Dribble: dash would be collision, but we can't kick the ball?"));
            Mem->LogAction2(70, "dribble: dash would collide, but can't kick!");
        }
    }

    // SMURF: need to plan for end point
    /* it might be a good idea to plan for the point where we want to stop
       dribbling so we don't overshoot and/or lose the ball when we stop */

    Mem->LogAction2(80, "dribble: doing a dash");
    DebugDrib(printf("Dribble: doing a basic dash\n"));
    com.time = Mem->CurrentTime;
    com.type = CMD_dash;
    com.power = max_pow;
    return com;
}

/* returns a kick that would put that would (assuming we dash next cycle)
   leave us with the ball at the correct angle after the next cycle */
TurnKickCommand GetDribbleKickCom(Vector vMyPred, AngleDeg drib_ang)
{
    Vector vBallTarg;

    Mem->LogAction2(90, "dribble: getting the kick command");

    if (!Mem->BallKickable())
        my_error("GetDribbleKickCom: Ball isn't kickable\n");

    vBallTarg = vMyPred +
                Polar2Vector(Mem->CP_dribble_ball_dist, drib_ang + Mem->MyBodyAng());

    vBallTarg = Mem->Global2RelativeToMyBody(vBallTarg);
    Vector vBallTraj = vBallTarg - Mem->BallRelativeToBodyPosition();
    /*vBallTraj = vBallTraj.rotate(-Mem->MyAng());*/

    return dokick(vBallTraj.dir(), vBallTraj.mod() / (1 + Mem->SP_ball_decay));
}

/* we may need to do a turnkick if the ball gets stuck behind us for some
   reason */
TurnKickCommand GetDribbleTurnKickCom(AngleDeg drib_ang)
{
    TurnKickCommand com;

    Mem->LogAction2(90, "dribble: getting the turnkick command");

    if (!Mem->BallKickable())
        my_error("GetDribbleTurnKickCom: Ball isn't kickable\n");

    DebugDrib(printf("Dribble: Having to do turnball\n"));
    TurnDir rot = RotToAvoidOpponent(drib_ang);
    DebugDrib(printf("Dribble: rot dir is %d\n", rot));
    KickToRes res = turnball_kick(drib_ang, rot, FALSE, &com);
    if (res == KT_LostBall)
        my_error("Dribble: turning lost the ball");
    else if (res == KT_Success || res == KT_DidNothing)
    {
        /* if we don;t know where we want to turn, just turn to in front of us
       so that our next kick won;t be a collision */
        res = turnball_kick(0, TURN_CLOSEST, FALSE, &com);
        if (res == KT_Success || res == KT_DidNothing)
            my_error("Dribble: TurnKick finished in dribble, so we shouldn't have doen it");
    }
    return com;
}

/* figures out a kick to do to keep dribbling
   (see the dribble_kick comment).
   Checks for collisions and calls turnkick if applicable */
/* the DribbleMode determines what to do when the kick we want to do is a
   collision. In strict mode, we always stop and turnball.
   In lazy mode, we try to get halfway to the desired angle and if that doesn't
   work, then we turnball */
/* lazy mode has been tested and works, but Strict mode has not. It seems to
   lead to many lost balls */
TurnKickCommand dribble_kick(Vector vEndPos, float max_pow,
                             AngleDeg drib_ang, DribbleMode mode)
{
    DebugDrib(cout << "in dribble_kick" << endl);
    if (!Mem->BallKickable())
        my_error("dribble_kick: Ball isn't kickable\n");

    TurnKickCommand com;
    /* figure out the point to kick the ball to */
    /* the 1 at the end is an idle cycle */
    Vector vMyPred = Mem->MyPredictedPosition(2, max_pow, 1);

    if ((vMyPred - Mem->MyPos()).mod() > (vEndPos - Mem->MyPos()).mod())
    {
        /* if we are going to go far enough to get to our endpos, we may not want
           to dash as far */
        vMyPred = vEndPos;
    }

    com = GetDribbleKickCom(vMyPred, drib_ang);

    if (Mem->WillKickBeCollision(com.power, com.angle, Mem->CP_collision_buffer))
    {
        if (mode == DM_Strict && fabs(Mem->BallAngleFromBody() - drib_ang) > 20)
        {
            //	(mode == DM_Lazy && signf(Mem->BallAngle()) == signf(drib_ang))) {
            /* we'll turnball instead */
            DebugDrib(printf("In strict mode, doing turnball\n"));
            Mem->LogAction2(80, "dribble: strict mode, turnballing");
            com = GetDribbleTurnKickCom(drib_ang);
        }
        else if (mode == DM_Lazy)
        {
            // SMURF: should we average, or send to 0???
            // com = GetDribbleKickCom(vMyPred, 0);
            DebugDrib(printf("Dribble: In lazy mode, trying to interpolate changeover\n"));
            com = GetDribbleKickCom(vMyPred,
                                    (Mem->BallAngleFromBody() + drib_ang) / 2);
            if (Mem->WillKickBeCollision(com.power, com.angle, Mem->CP_collision_buffer))
            {
                /* still didn't work, use turnball! */
                DebugDrib(printf("Dribble: In lazy mode, didn't work\n"));
                Mem->LogAction2(80, "dribble: lazy mode, but still have to turnball");
                com = GetDribbleTurnKickCom(drib_ang);
            }
            else
            {
                Mem->LogAction2(80, "dribble: lazy mode, averaging kick");
            }
        }
        else
            my_error("Dribble: mode is bad!");
    }
    else
    {
        DebugDrib(printf("Dribble: doing a kick\n"));
        Mem->LogAction2(80, "dribble: doing a normal kick");
    }

    return com;
}

DribbleRes DribbleTo(Vector vEndPos, float max_dash_pow, float buffer,
                     AngleDeg drib_ang, DribbleMode mode,
                     Bool IsDodge, Vector DodgePoint)
{
    AngleDeg targAng;
    AngleDeg targAngErr = Mem->CP_max_go_to_point_angle_err;

    if (!Mem->BallPositionValid())
    {
        my_error("must know where ball is to dribble");
        return DR_Error;
    }

    /* first check velocity confidence */
    if (!Mem->BallVelocityValid())
    {
        my_error("Shouldn't call dribble with ball velocity invalid");
        return DR_Error;
        /*Mem->LogAction2(60, "Dribble: ball velocity not valid");
        face_only_neck_to_ball();
        stop_ball();
        return DR_Going; */
    }

    DebugDrib(cout << "MyPos: " << Mem->MyPos() << "\tMyAng: " << Mem->MyAng() << endl);
    DebugDrib(cout << "BallPos: " << Mem->BallAbsolutePosition() << endl);
    DebugDrib(cout << "BallDist: " << Mem->BallDistance() << endl);
    DebugDrib(cout << "BallVel: " << Mem->BallAbsoluteVelocity() << endl);
    DebugDrib(cout << "vEndPos: " << vEndPos << "\tdrib_ang: " << drib_ang << endl);

    Mem->LogAction6(60, "DribbleTo: pow %.1f to (%.1f, %.1f) at ang %.2f",
                    max_dash_pow, vEndPos.x, vEndPos.y, drib_ang);

    if (mode == DM_None)
        my_error("Can't pass DM_None as a mode to DribbleTo");

    if ((Mem->MyPos() - vEndPos).mod() < buffer)
    {
        DebugDrib(printf("Got to target!\n"));
        Mem->LogAction2(80, "dribble: I got to target!");
        return DR_GotThere;
    }

    if (IsDodge)
    {
        /* there is something we need to dodge */
        /* we'll dodge whatever side of us the ball was supposed to go */
        Vector disp = vEndPos - Mem->MyPos();
        // DebugDrib(cout << "Disp1: " << di	sp << endl);
        disp *= Mem->CP_dodge_distance_buffer / disp.mod();
        // DebugDrib(cout << "Disp2: " <<	 disp << endl);
        disp = disp.rotate(signf(drib_ang) * 90);
        DebugDrib(cout << "Disp: " << disp << endl);
        if (Mem->DistanceTo(DodgePoint) < Mem->CP_dribble_dodge_close_dist)
        {
            DebugDrib(printf("Dribble: sharp dodge\n"));
            vEndPos = Mem->MyPos() + disp;
        }
        else
        {
            DebugDrib(printf("Dribble: regular dodge\n"));
            vEndPos = (DodgePoint + disp - Mem->MyPos()) * 2 + Mem->MyPos();
        }
        DebugDrib(cout << "We're dodging a player: disp: " << disp
                       << "\tvEndPos: " << vEndPos << endl);
        Mem->LogAction2(80, "dribble: dodging a player");
        targAngErr = Mem->CP_dribble_dodge_angle_err;
    }

    targAng = Mem->AngleToFromBody(vEndPos);
    DebugDrib(printf("targAng: %f\n", targAng));
    TurnKickCommand com;

    /* if the mode is strict and then ball is way off,
       then go ahead and turnball */
    if (mode == DM_Strict &&
        fabs(Mem->BallAngleFromBody() - (targAng + drib_ang)) >
            Mem->CP_dribble_exp_angle_buffer &&
        Mem->BallKickable())
    {
        com = GetDribbleTurnKickCom(drib_ang + targAng);
        if (DoTurnKickCommand(com))
            return DR_Going;
        else
            return DR_Error;
    }

    /* first turn the right way */
    if (fabs(targAng) > targAngErr)
    {
        DebugDrib(printf("Want to turn\n"));
        if (!Mem->BallWillBeKickable(1, 0, Mem->CP_kickable_buffer) &&
            !Mem->BallWillBeKickable(2, 0, 0))
            if (Mem->BallKickable())
            {
                DebugDrib(printf("Ball will not be kickable, but it is now, so we're kicking\n"));
                Mem->LogAction2(80, "dribble: stopping the ball");
                com = dokick(0, 0);
            }
            else
            {
                DebugDrib(printf("Ball will not be kickable, trying to catch up\n"));
                Mem->LogAction2(70, "dribble: trying to catch up to ball");
                com = dribble_dash(vEndPos, max_dash_pow, drib_ang, mode);
                if (com.time != Mem->CurrentTime)
                    return DR_LostBall;
            }
        else
        {
            DebugDrib(printf("Dribble: turning to the right angle %f %f\n", targAng, Mem->MyAng()));
            Mem->LogAction2(80, "dribble: turning the right way");
            com.time = Mem->CurrentTime;
            com.type = CMD_turn;
            com.angle = targAng;
            com.turn_neck = FALSE;
        }
        if (DoTurnKickCommand(com))
            return DR_Going;
        else
            return DR_Error;
    }

    if (Mem->BallKickable())
    {
        AngleDeg expBallAngle;
        expBallAngle = (Mem->BallPredictedPosition() -
                        Mem->MyPredictedPosition(1, max_dash_pow))
                           .dir();
        expBallAngle -= Mem->MyBodyAng();
        NormalizeAngleDeg(&expBallAngle);
        DebugDrib(printf("Expected ball angle: %f\n", expBallAngle));
        if (fabs(expBallAngle - drib_ang) > Mem->CP_dribble_exp_angle_buffer ||
            fabs(expBallAngle) > 90)
        {
            /* kick if the ball will end up behind us */
            DebugDrib(printf("The ball will be too far off it's target, so kicking\n"));
            Mem->LogAction2(70, "dribble: kicking because ball will be behind us");
            com = dribble_kick(vEndPos, max_dash_pow, drib_ang, mode);
        }
        else
        {
            /* otherwise dash */
            DebugDrib(printf("Going to try to do a dash\n"));
            Mem->LogAction2(70, "dribble: trying to dash");
            com = dribble_dash(vEndPos, max_dash_pow, drib_ang, mode);
        }
    }
    else
    {
        DebugDrib(printf("Ball is not kickable, better dash\n"));
        Mem->LogAction2(70, "dribble: trying to dash because ball not kickable");
        com = dribble_dash(vEndPos, max_dash_pow, drib_ang, mode);
    }

    if (com.time != Mem->CurrentTime)
        return DR_LostBall;

    if (DoTurnKickCommand(com))
        return DR_Going;
    else
        return DR_Error;
}

DribbleRes SimpleDribbleTo(Vector vEndPos, float max_dash_pow, float buffer)
{
    return DribbleTo(vEndPos, max_dash_pow, buffer, Mem->CP_dribble_angle_norm,
                     DM_Lazy, FALSE, Vector(0, 0));
}

AngleDeg GetNormalizeDribAng(AngleDeg ang)
{
    NormalizeAngleDeg(&ang);
    if (ang <= -90)
        ang = -90;
    if (ang >= 90)
        ang = 90;
    return ang;
}

AngleDeg GetNormalizeDribAngWithBuffer(AngleDeg ang)
{
    NormalizeAngleDeg(&ang);
    if (fabs(ang) > 180 - Mem->CP_dribble_angle_ignore_buffer)
        /* we want the ball right in back, so use the side the ball is on */
        ang = signf(Mem->BallAngleFromBody()) * 90;
    else
        ang = GetNormalizeDribAng(ang);
    return ang;
}

/* -*- Mode: C++ -*- */

/* behave.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

#define DebugMark(x)

#ifdef DEBUG_OUTPUT
#define DebugBlock(x) x
#else
#define DebugBlock(x)
#endif

//#define MAKE_DEMO
//#define TEST_PETER
//#define TEST_PAT
//#define TEST_PAT2

void behave()
{

    /* Put any one of these test functions in -- but only one at a time */

    /* Start one player -- move the ball around the field */
    // face_ball();

    /* Start one player, press kickoff -- it will kick to the goal */
    test_go_to_ball();

    /* Start one player, press kickoff -- it will dribble to the goal */
    // test_dribble_to_goal();

    /* Start 2 players on different teams, press 'drop ball' near one player */
    // test_1v1();

    /* Start 2 players on different teams, press 'drop ball' near one player */
    // test_volley();

    /* These are to print out the positions of objects */
    // test_print_ball();
    // test_print_positions();

    /* Start 3 players on each team and press 'drop ball' in the middle of the field */
    // Rectangle *rect = new Rectangle(Vector(0,0),Vector(30,30));
    // test_random_movement_in_rectangle(rect);
    // delete rect;

    /* make sure that the option save_action_log_level is >= 200
       The player will generate a file
       Logfiles/<team_name><num>-<side>-actions.log which records state
       information and action choices. This logfile is pure text. It can
       be viewed more easily using the modified logplayer which is
       available with the CMUnited99 support code. This is part of a
       system we call Layered Extropsection. Please see the web page
       indicated above for more details */
    // test_log_action();

    /* Look in test.[Ch] for more interesting beahviors to play with */
}

/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/

ActionQueueRes scan_field_with_body()
{
    Mem->LogAction3(40, "scan_field_with_body (%d)", Mem->TimeToTurnForScan());
    if (Mem->TimeToTurnForScan())
    {
        turn(Mem->MyViewAngle() * 2 - Mem->CP_scan_overlap_angle);
        return AQ_ActionQueued;
    }
    else
        return AQ_ActionNotQueued;
}

/*****************************************************************************************/

void turn_neck_to_relative_angle(AngleDeg ang)
{
    turn_neck(GetNormalizeAngleDeg(ang - Mem->MyNeckRelAng()));
}

/*****************************************************************************************/

void scan_field_with_neck()
{
    Mem->LogAction3(40, "scan_field_with_neck (%d)", Mem->TimeToTurnForScan());
    if (Mem->TimeToTurnForScan())
    {
        if (Mem->MyNeckRelAng() >= Mem->SP_max_neck_angle - 1) /* take into account reporting error */
            turn_neck_to_relative_angle(Mem->SP_min_neck_angle);
        else
            turn_neck(Mem->LimitTurnNeckAngle(Mem->MyViewAngle() * 2 - Mem->CP_scan_overlap_angle));
    }
}

/*****************************************************************************************/

ActionQueueRes face_only_body_to_point(Vector point)
{
    /* don't turn neck */
    Mem->LogAction4(30, "facing only body to point (%.1f %.1f)", point.x, point.y);

    /* shouldn't actually have queued actions at this point */
    AngleDeg target_rel_ang = Mem->PredictedPointRelAngFromBodyWithQueuedActions(point);

    if (fabs(target_rel_ang) < 1)
    {
        Mem->LogAction2(40, "Already close enough");
        return AQ_ActionNotQueued;
    }
    turn(target_rel_ang);
    return AQ_ActionQueued;
}

/*****************************************************************************************/

void face_only_neck_to_point(Vector point)
{
    /* don't turn body */
    Mem->LogAction4(30, "facing only neck to point (%.1f %.1f)", point.x, point.y);

    AngleDeg target_rel_ang = Mem->PredictedPointRelAngFromBodyWithQueuedActions(point);

    if (fabs(GetNormalizeAngleDeg(Mem->MyNeckRelAng() - target_rel_ang)) < 1)
    {
        Mem->LogAction2(40, "Already close enough");
        return;
    }

    if (Mem->CanSeeAngleFromBodyWithNeck(target_rel_ang))
    {
        turn_neck(Mem->LimitTurnNeckAngle(target_rel_ang - Mem->MyNeckRelAng()));
    }
    else
        Mem->LogAction5(30, "can't face point (%.1f %.1f) with only neck (%.1f)",
                        point.x, point.y, target_rel_ang);
}

/*****************************************************************************************/

ActionQueueRes face_neck_to_point(Vector point)
{
    /* face_neck can turn body if needed */
    Mem->LogAction4(30, "facing neck to point (%.1f %.1f)", point.x, point.y);

    AngleDeg target_rel_ang = Mem->PredictedPointRelAngFromBodyWithQueuedActions(point);

    if (fabs(GetNormalizeAngleDeg(Mem->MyNeckRelAng() - target_rel_ang)) < 1)
    {
        Mem->LogAction2(40, "Already close enough");
        return AQ_ActionNotQueued;
    }

    if (Mem->CanFaceAngleFromBodyWithNeck(target_rel_ang))
    {
        Mem->LogAction2(35, "can face with neck");
        turn_neck_to_relative_angle(target_rel_ang);
        return AQ_ActionNotQueued;
    }

    /* If can't do it with just neck, turn body as much as needed to face directly */
    AngleDeg max_turn = Mem->MaxEffectiveTurn();
    if (fabs(target_rel_ang) < max_turn)
    {
        Mem->LogAction2(35, "can't face with neck, can with body");
        turn(target_rel_ang);
        turn_neck_to_relative_angle(0);
        return AQ_ActionQueued;
    }

    Mem->LogAction2(35, "can't face with neck or body alone, turning both");
    turn(target_rel_ang);
    target_rel_ang -= max_turn; /* The neck target_ang */

    if (target_rel_ang < Mem->SP_min_neck_angle)
    {
        Mem->LogAction2(40, "couldn't face all the way");
        turn_neck_to_relative_angle(Mem->SP_min_neck_angle);
    }
    else if (target_rel_ang > Mem->SP_max_neck_angle)
    {
        Mem->LogAction2(40, "couldn't face all the way");
        turn_neck_to_relative_angle(Mem->SP_max_neck_angle);
    }
    else
        turn_neck_to_relative_angle(target_rel_ang);

    return AQ_ActionQueued;
}

/*****************************************************************************************/

ActionQueueRes face_neck_and_body_to_point(Vector point)
{
    /* face_neck_and_body will turn both as much as possible to the point */
    Mem->LogAction4(30, "facing neck and body to point (%.1f %.1f)", point.x, point.y);

    AngleDeg max_turn = Mem->MaxEffectiveTurn();
    AngleDeg target_rel_ang = Mem->PredictedPointRelAngFromBodyWithQueuedActions(point);

    if (fabs(GetNormalizeAngleDeg(Mem->MyNeckRelAng() - target_rel_ang)) < 1 && fabs(target_rel_ang) < 1)
    {
        Mem->LogAction2(40, "Already close enough");
        return AQ_ActionNotQueued;
    }

    if (fabs(target_rel_ang) < max_turn)
    {
        Mem->LogAction2(35, "Can get both neck and body there");
        /* Can get both neck and body there */
        face_only_body_to_point(point);
        turn_neck_to_relative_angle(0);
        return AQ_ActionQueued;
    }

    /* Turn body as much as possible and try to get neck there */
    return face_neck_to_point(point);
}

/*****************************************************************************************/

ActionQueueRes face_only_body_to_player(char side, Unum num)
{
    if (Mem->PlayerPositionValid(side, num))
    {
        return face_only_body_to_point(Mem->PlayerAbsolutePosition(side, num));
    }
    else
        return scan_field_with_body();
}

/*****************************************************************************************/

void face_only_neck_to_player(char side, Unum num)
{
    if (Mem->PlayerPositionValid(side, num))
    {
        face_only_neck_to_point(Mem->PlayerAbsolutePosition(side, num));
    }
    else
        scan_field_with_neck();
}

/*****************************************************************************************/

ActionQueueRes face_neck_to_player(char side, Unum num)
{
    if (Mem->PlayerPositionValid(side, num))
    {
        return face_neck_to_point(Mem->PlayerAbsolutePosition(side, num));
    }
    else
        return scan_field_with_body();
}

/*****************************************************************************************/

ActionQueueRes face_neck_and_body_to_player(char side, Unum num)
{
    if (Mem->PlayerPositionValid(side, num))
    {
        return face_neck_and_body_to_point(Mem->PlayerAbsolutePosition(side, num));
    }
    else
        return scan_field_with_body();
}

/*****************************************************************************************/

ActionQueueRes face_only_body_to_opponent(Unum opponent)
{
    Mem->LogAction3(30, "facing only body to opponent %d", opponent);
    return face_only_body_to_player(Mem->TheirSide, opponent);
}

/*****************************************************************************************/

void face_only_neck_to_opponent(Unum opponent)
{
    Mem->LogAction3(30, "facing only neck to opponent %d", opponent);
    face_only_neck_to_player(Mem->TheirSide, opponent);
}

/*****************************************************************************************/

ActionQueueRes face_neck_to_opponent(Unum opponent)
{
    Mem->LogAction3(30, "facing neck to opponent %d", opponent);
    return face_neck_to_player(Mem->TheirSide, opponent);
}

/*****************************************************************************************/

ActionQueueRes face_neck_and_body_to_opponent(Unum opponent)
{
    Mem->LogAction3(30, "facing neck and body to opponent %d", opponent);
    return face_neck_and_body_to_player(Mem->TheirSide, opponent);
}

/*****************************************************************************************/

ActionQueueRes face_only_body_to_teammate(Unum teammate)
{
    Mem->LogAction3(30, "facing only body to teammate %d", teammate);
    return face_only_body_to_player(Mem->MySide, teammate);
}

/*****************************************************************************************/

void face_only_neck_to_teammate(Unum teammate)
{
    Mem->LogAction3(30, "facing only neck to teammate %d", teammate);
    face_only_neck_to_player(Mem->MySide, teammate);
}

/*****************************************************************************************/

ActionQueueRes face_neck_to_teammate(Unum teammate)
{
    Mem->LogAction3(30, "facing neck to teammate %d", teammate);
    return face_neck_to_player(Mem->MySide, teammate);
}

/*****************************************************************************************/

ActionQueueRes face_neck_and_body_to_teammate(Unum teammate)
{
    Mem->LogAction3(30, "facing neck and body to teammate %d", teammate);
    return face_neck_and_body_to_player(Mem->MySide, teammate);
}

/*****************************************************************************************/

ActionQueueRes face_only_body_to_ball()
{
    Mem->LogAction2(30, "facing body to ball");
    if (Mem->BallPositionValid())
    {
        return face_only_body_to_point(Mem->BallPredictedPositionWithQueuedActions());
    }
    else
        return scan_field_with_body();
}

/*****************************************************************************************/

void face_only_neck_to_ball()
{
    Mem->LogAction2(30, "facing only neck to ball");
    if (Mem->BallPositionValid())
    {
        face_only_neck_to_point(Mem->BallPredictedPositionWithQueuedActions());
    }
    else
        scan_field_with_neck();
}

/*****************************************************************************************/

ActionQueueRes face_neck_to_ball()
{
    Mem->LogAction2(30, "facing neck to ball");
    if (Mem->BallPositionValid())
    {
        return face_neck_to_point(Mem->BallPredictedPositionWithQueuedActions());
    }
    else
        return scan_field_with_body();
}

/*****************************************************************************************/

ActionQueueRes face_neck_and_body_to_ball()
{
    Mem->LogAction2(30, "facing neck and body to ball");
    if (Mem->BallPositionValid())
    {
        return face_neck_and_body_to_point(Mem->BallPredictedPositionWithQueuedActions());
    }
    else
        return scan_field_with_body();
}

/*****************************************************************************************/
/* if the arg is DT_all, we always dodge, otherwise we only dodge if they don't have ball */
void get_ball()
{
    if (!Mem->MyConf() || !Mem->BallPositionValid())
        my_error("not enough info to get ball");

    DodgeType dodge = DT_only_with_ball;

    if (!Mem->BallMoving())
    {
        Mem->LogAction2(30, "get_ball: ball not moving, going to it's pos");
        if (go_to_point(Mem->BallAbsolutePosition(), 0, 100, dodge) == AQ_ActionNotQueued)
        {
            my_error("already there???");
            face_neck_and_body_to_ball();
        }
        face_only_neck_to_ball();
    }
    else
    {
        if (!Mem->MyInterceptionAble())
        {
            Mem->LogAction2(30, "get_ball: going to the moving ball, but can't?");
            my_error("Can't get to the ball");
            face_neck_and_body_to_ball();
        }
        else if (Mem->MyInterceptionNumberCycles() == 1)
        {
            /* we're just one dash away, so just do it */
            Mem->LogAction2(30, "get_ball: going to the moving ball, just dashing 1 cycle");
            dash(Mem->CorrectDashPowerForStamina(Mem->MyInterceptionDashPower()));
            face_only_neck_to_ball();
        }
        else if (go_to_point(Mem->MyInterceptionPoint(), 0,
                             Mem->MyInterceptionDashPower(), dodge) == AQ_ActionNotQueued)
        {
            Mem->LogAction2(30, "get_ball: going to the moving ball, but already there?");
            my_error("already there (moving) ???");
            face_neck_and_body_to_ball();
        }
        else
        {
            Mem->LogAction4(30, "get_ball: going to the moving ball (%d) pow %.2f",
                            Mem->MyInterceptionNumberCycles(), Mem->MyInterceptionDashPower());
            face_only_neck_to_ball();
        }
    }
}

/*****************************************************************************************/

void stop_ball()
{
    if (!Mem->BallPositionValid())
        my_error("Need to know where ball is to stop it");

    if (Mem->BallVelocityValid())
    {
        Mem->LogAction2(30, "stop_ball: velocity valid");
        DoTurnKickCommand(dokick(0, 0));
    }
    else
    {
        Mem->LogAction2(30, "stop_ball: velocity not valid");
        // DebugDrib(printf("Stop kick; don't know velocity, so goign to kick to us"));
        kick(Mem->CP_stop_ball_power, GetNormalizeAngleDeg(Mem->BallAngleFromBody() + 180));
    }
}

/*****************************************************************************************/

void hold_ball()
{
    if (!Mem->BallPositionValid())
        my_error("Need to know where ball is to hold it");

    Unum opponent = Mem->ClosestOpponent();

    Mem->LogAction3(30, "hold_ball: closest opponent == %d", opponent);

    /* make sure that we keep control of the ball
       we need to make sure that ball's velocity is valid,
       and then scan_field only if we will keep control of ball */

    if (opponent == Unum_Unknown)
    {
        if (!Mem->BallVelocityValid())
        {
            Mem->LogAction2(40, "hold_ball: velocity not valid");
            face_neck_to_ball();
        }
        else if (Mem->TimeToTurnForScan() &&
                 (!Mem->BallMoving() ||
                  Mem->BallWillBeKickable(1, 0, Mem->CP_holdball_kickable_buffer)))
        {
            Mem->LogAction2(40, "hold_ball: ball not moving or will be kickable");
            scan_field_with_body();
        }
        else if (Mem->BallMoving() ||
                 !Mem->BallWillBeKickable(1, 0, Mem->CP_holdball_kickable_buffer))
        {
            Mem->LogAction2(40, "hold_ball: ball moving or won't be kickable");
            stop_ball();
        }
        else
        {
            Mem->LogAction2(40, "hold_ball: doing nothing");
            /* do nothing */
        }
    }
    else
    {
        /* closest opponent known */
        AngleDeg ang = GetNormalizeAngleDeg(Mem->OpponentAngleFromBody(opponent) + 180);

        Vector targ_pos = Mem->RelativeToMyBody2Global(Polar2Vector(Mem->CP_opt_ctrl_dist, ang));
        if (!Mem->FieldRectangle.IsWithin(targ_pos))
        {
            /* Adjust the targ_pos and ang to be in bounds- ignore the corners! */
            Line lSide = Mem->FieldRectangle.nearestEdgeLine(Mem->MyPos());
            Vector sol1, sol2;
            int num_sol = LineCircleIntersect(lSide, Mem->CP_opt_ctrl_dist, Mem->MyPos(),
                                              &sol1, &sol2);
            if (num_sol == 0)
            {
                my_error("hold_ball: why didn't LineCircleIntersect work?");
            }
            else if (num_sol == 1)
            {
                Mem->LogAction6(40, "hold_ball: avoiding the sideline 1; old: (%.1f, %.1f)  new: (%.1f, %.1f)",
                                targ_pos.x, targ_pos.y, sol1.x, sol1.y);
                ang = Mem->AngleToFromBody(sol1);
            }
            else if (num_sol == 2)
            {
                if (targ_pos.dist2(sol1) < targ_pos.dist2(sol2))
                {
                    Mem->LogAction6(40, "hold_ball: avoiding the sideline 2; old: (%.1f, %.1f)  new: (%.1f, %.1f)",
                                    targ_pos.x, targ_pos.y, sol1.x, sol1.y);
                    ang = Mem->AngleToFromBody(sol1);
                }
                else
                {
                    Mem->LogAction6(40, "hold_ball: avoiding the sideline 3; old: (%.1f, %.1f)  new: (%.1f, %.1f)",
                                    targ_pos.x, targ_pos.y, sol2.x, sol2.y);
                    ang = Mem->AngleToFromBody(sol2);
                }
            }
        }

        if (!Mem->BallWillBeKickable(1, 0, Mem->CP_holdball_kickable_buffer))
        {
            Mem->LogAction2(40, "hold_ball: turning ball from opponent -- won't be kickable");
            if (TurnballTo(ang) == KT_DidNothing)
                my_error("hold_ball: why didn't turnball do something");
        }
        else if (Mem->OpponentPositionValid(opponent) < .9 &&
                 !Mem->InViewAngle(Mem->OpponentAngleFromNeck(opponent)))
        {
            Mem->LogAction2(40, "hold_ball: looking at opponent");
            face_neck_to_opponent(opponent);
        }
        else if (Mem->TimeToTurnForScan() && Mem->EstimatedCyclesToSteal(opponent) > Mem->CP_time_for_full_rotation / 3)
        {
            /* do we want to scan when the ball may be near the opp? */
            /* yes, but only if the opponent isn't about to get the ball */
            Mem->LogAction2(40, "hold_ball: looking around");
            scan_field_with_body();
        }
        else
        {
            Mem->LogAction2(40, "hold_ball: turning ball from opponent");
            TurnballTo(ang);
        }
    }
}

/*****************************************************************************************/

void pass_ball(Unum teammate, float target_vel_at_dest)
{
    Mem->LogAction3(30, "passing to %d", teammate);

    if (teammate == Unum_Unknown)
        my_error("Need to pass to a teammate");
    if (!Mem->TeammatePositionValid(teammate))
        my_error("can't pass to invalid teammmate");

    if (Mem->TeammatePositionValid(teammate) < .9)
    {
        Mem->LogAction2(40, "pass_ball: looking for teammate");
        if (face_neck_to_teammate(teammate) == AQ_ActionNotQueued)
            hold_ball();
        return;
    }

    Vector target = Mem->TeammateAbsolutePosition(teammate);
    AngleDeg target_angle = Mem->AngleToFromBody(target);
    TurnDir rotation = KickRotationDirection(target_angle);
    float target_vel = Mem->VelAtPt2VelAtFoot(target, target_vel_at_dest);
    KickMode mode = Mem->BestKickMode(target_angle);
    if (mode == KM_HardestKick)
        mode = KM_Hard;

    Mem->team_passer = Mem->MyNumber;
    Mem->team_receiver = teammate;
    Mem->team_pass_time = Mem->CurrentTime;

    Mem->LogAction2(40, "pass_ball: starting pass");
    kick_ball(target, mode, target_vel, rotation);
}

/*****************************************************************************************/

/* extend this angle to the sideline and call the other kick_ball variant */
/* target_angle is relative to body */
void kick_ball(AngleDeg target_angle, KickMode mode, float target_vel, TurnDir rotation)
{
    Vector target_pt =
        Mem->FieldRectangle.RayIntersection(Ray(Mem->MyPos(), target_angle + Mem->MyBodyAng()));
    Mem->LogAction5(40, "starting kick to angle %.1f, translated to point (%.1f, %.1f)",
                    target_angle, target_pt.x, target_pt.y);
    kick_ball(target_pt, mode, target_vel, rotation);
    return;
}

void kick_ball(Vector point, KickMode mode, float target_vel, TurnDir rotation)
{
    if (rotation == TURN_NONE)
        rotation = KickRotationDirection(Mem->AngleToFromBody(point));

    /* look to see if a dash will help */
    if (mode == KM_Hard && Mem->WillDashHelpKick(point, Mem->SP_max_power))
    {
        mode = KM_Moderate;
        float target_angle =
            (point - Mem->MyPredictedPosition(1, Mem->SP_max_power)).dir() - Mem->MyBodyAng();
        Mem->LogAction3(40, "kick_ball: dashing will help kick to angle: %1.f", target_angle);
        dash(Mem->SP_max_power);
        Mem->StartKick(target_angle, mode, target_vel, rotation);
    }
    else
    {
        if (mode == KM_Hard)
            mode = KM_Moderate;
        float target_angle = Mem->AngleToFromBody(point);
        Mem->LogAction3(40, "kick_ball: starting kick to angle %.1f", target_angle);
        Mem->StartKick(target_angle, mode, target_vel, rotation);
        smart_kick_hard(target_angle, mode, target_vel, rotation);
    }

    return;
}

void kick_ball(AngleDeg target_angle, KickMode mode, TurnDir rotation)
{
    kick_ball(target_angle, mode, 2 * Mem->SP_ball_speed_max, rotation);
}

void kick_ball(Vector point, KickMode mode, TurnDir rotation)
{
    kick_ball(point, mode, 2 * Mem->SP_ball_speed_max, rotation);
}

/*****************************************************************************************/

ActionQueueRes go_to_point(Vector p, float buffer, float dash_power, DodgeType dodge)
{
    Mem->LogAction5(30, "go_to_point %d (%.1f %.1f)", dodge, p.x, p.y);
    if (!Mem->MyConf())
        my_error("Can't go to a point if not localized");

    if (Mem->DistanceTo(p) < buffer)
    {
        if (Mem->SP_use_offside && fabs(Mem->MyX() - Mem->my_offside_line) < 5)
        { /* hack */
            Unum opp = Mem->FurthestForwardOpponent();
            if (opp != Unum_Unknown && Mem->OpponentPositionValid(opp) < .9)
            { /* hack */
                Mem->LogAction2(40, "go_to_point: looking for offsides line");
                return face_neck_to_opponent(opp);
            }
        }
        Mem->LogAction2(40, "go_to_point: already at the point");
        return AQ_ActionNotQueued;
    }

    if (Mem->PlayMode == PM_Their_Goal_Kick && Mem->MyPos() != p)
    {
        /* if ( Mem->TheirPenaltyArea.IsWithin(p) ){
          my_error("Can't go into their penalty area on a goal kick!"); */
        Line l = LineFromTwoPoints(Mem->MyPos(), p);
        Vector intersection = AdjustPtToRectOnLine(Mem->MyPos(), Mem->TheirPenaltyArea, l);
        if (intersection != Mem->MyPos() && l.InBetween(intersection, Mem->MyPos(), p))
        {
            /* Need to go around the rectangle */
            Mem->LogAction2(40, "go_to_point: moving around penalty area");
            Vector target;
            if (Mem->MyX() < Mem->TheirPenaltyArea.LeftX())
                target = Vector(Mem->TheirPenaltyArea.LeftX() - 3, p.y);
            else if (Mem->MyY() > 0)
                target = Vector(Mem->TheirPenaltyArea.LeftX() - 3, Mem->TheirPenaltyArea.BottomY() + 3);
            else
                target = Vector(Mem->TheirPenaltyArea.LeftX() - 3, Mem->TheirPenaltyArea.TopY() - 3);
            go_to_point(target, 0, dash_power, dodge);
            return AQ_ActionQueued;
        }
    }

    if (Mem->PlayMode != PM_Play_On && Mem->TeamInPossession() == Mem->TheirSide &&
        !Mem->OffsidePosition(p, Mem->MySide) &&
        // p.dist(Mem->BallAbsolutePosition()) > Mem->SP_free_kick_buffer &&
        Mem->InOffsidePosition() && Mem->BallDistance() < Mem->SP_free_kick_buffer + 1)
    {
        Mem->LogAction2(40, "go_to_point: moving around free_kick area");
        if (Mem->BallY() > Mem->MyY())
            go_to_point(Vector(Mem->MyX(), Mem->BallY() - (Mem->SP_free_kick_buffer + 1)));
        else
            go_to_point(Vector(Mem->MyX(), Mem->BallY() + (Mem->SP_free_kick_buffer + 1)));
        return AQ_ActionQueued;
    }

    float target_ang = GetNormalizeAngleDeg((p - Mem->MyPredictedPosition()).dir() - Mem->MyBodyAng());
    float target_dist = Mem->DistanceTo(p);

    if (dodge != DT_none)
    { /* dodge players */
        PlayerObject *player;
        float dodge_dist = Min(Mem->CP_dodge_distance_buffer, target_dist);
        AngleDeg dodge_ang = Mem->CP_dodge_angle_buffer;
        if ((player = Mem->GetPlayerWithin(dodge_dist, dodge_ang, 0, target_ang - dodge_ang)) != NULL &&
            (dodge != DT_unless_with_ball ||
             (Mem->BallPositionValid() &&
              player->get_abs_pos().dist(Mem->BallAbsolutePosition()) > Mem->SP_kickable_area)) &&
            (dodge != DT_only_with_ball ||
             (Mem->BallPositionValid() &&
              player->get_abs_pos().dist(Mem->BallAbsolutePosition()) <= Mem->SP_kickable_area)))
        {
            Mem->LogAction2(40, "go_to_point: dodging right");
            /*if ( Mem->NumPlayersWithin( dodge_dist, 2*dodge_ang) ){*/
            /* Target at dist player_size, so no players will be within in the next iteration ==> dash */
            Vector new_target = Mem->BodyPolar2Gpos(Mem->SP_player_size, player->get_ang_from_body() + Mem->CP_dodge_angle);
            if (new_target == p)
                my_error("Dodging isn't changing the point!");
            go_to_point(new_target, 0, Mem->CP_dodge_power, DT_none);
            /*}
            else{
          dash(Mem->CorrectDashPowerForStamina(Min(dash_power,Mem->CP_dodge_power)));
            }*/
            return AQ_ActionQueued;
        }
        if ((player = Mem->GetPlayerWithin(dodge_dist, dodge_ang, 0, target_ang + dodge_ang)) != NULL &&
            (dodge != DT_unless_with_ball ||
             (Mem->BallPositionValid() &&
              player->get_abs_pos().dist(Mem->BallAbsolutePosition()) > Mem->SP_kickable_area)) &&
            (dodge != DT_only_with_ball ||
             (Mem->BallPositionValid() &&
              player->get_abs_pos().dist(Mem->BallAbsolutePosition()) <= Mem->SP_kickable_area)))
        {
            Mem->LogAction2(40, "go_to_point: dodging left");
            /*if ( Mem->NumPlayersWithin( dodge_dist, 2*dodge_ang) ){*/
            /* Target at dist player_size, so no players will be within in the next iteration ==> dash */
            Vector new_target = Mem->BodyPolar2Gpos(Mem->SP_player_size, player->get_ang_from_body() - Mem->CP_dodge_angle);
            if (new_target == p)
                my_error("Dodging isn't changing the point!");
            go_to_point(new_target, 0, Mem->CP_dodge_power, DT_none);
            /*}
            else{
          dash(Mem->CorrectDashPowerForStamina(Min(dash_power,Mem->CP_dodge_power)));
            }*/
            return AQ_ActionQueued;
        }
    }

    if (fabs(target_ang) > Mem->CP_max_go_to_point_angle_err ||
        (Mem->PlayMode == PM_Their_Goal_Kick &&
         Mem->TheirPenaltyArea.IsWithin(Mem->MyPredictedPosition(1, dash_power))))
    {
        Mem->LogAction3(50, "go_to_point: turning %f", target_ang);
        turn(target_ang);
        return AQ_ActionQueued;
    }

    dash_power = Mem->CorrectDashPowerForStamina(dash_power);
    if (dash_power > 0)
    {
        Mem->LogAction3(50, "go_to_point: dashing %f", dash_power);
        dash(dash_power);
        return AQ_ActionQueued;
    }
    else
    {
        my_stamp;
        printf("recovering\n");
    }

    Mem->LogAction2(50, "go_to_point: doing nothing?");
    return AQ_ActionNotQueued;
}

/* -*- Mode: C++ -*- */

/* parse.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

void Parse_Sight(Time time, char *SightInfo);
void Parse_Sense(Time time, char *SenseInfo);
void Parse_Sound(Time time, char *SoundInfo);
void Parse_Referee_Sound(char *RefereeSound);
void Parse_Trainer_Sound(char *msg);
void Parse_My_Coach_Sound(Time time, char *msg);
void Parse_Their_Coach_Sound(Time time, char *msg);

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void Parse(char *SensoryInfo)
{
    SenseType sense_type;
    int time;
    switch (SensoryInfo[3])
    {
    case 'e':
        sense_type = See_Msg;
        break; /* see   */
    case 'n':
        sense_type = Sense_Msg;
        break; /* sense */
    case 'a':
        sense_type = Hear_Msg;
        break; /* hear  */
    default:
        my_error("Sent an illegal message");
        return;
    }

    time = get_int(&SensoryInfo); /* %d    */

    Time tm = Mem->update_time(time);

    switch (sense_type)
    {
    case See_Msg:
        if (!Mem->LastActionOpTime)
            break; /* Don't parse until I've started counting time   */
        if (tm == Mem->LastSightTime)
            break; /* Don't parse a second sight from the same cycle */
        if (Mem->NewSight == TRUE)
        {
            Mem->ClearSeenInfo();
            Mem->LogAction2(190, "Sight from last cycle lying around -- clearing it");
        }
        Parse_Sight(tm, SensoryInfo);
        Mem->LastSightInterval = tm - Mem->LastSightTime;
        Mem->LastSightTime = tm;
        Mem->NewSight = TRUE;
        break;
    case Sense_Msg:
        Parse_Sense(tm, SensoryInfo);
        Mem->LastSenseTime = tm;
        break;
    case Hear_Msg:
        if (!Mem->LastActionOpTime)
            break; /* Don't parse until I've started counting time */
        Parse_Sound(tm, SensoryInfo);
        Mem->LastSoundTime = tm;
        break;
    }

    Mem->LastSenseType = sense_type;
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void Parse_Sense(Time time, char *SenseInfo)
{
    get_word(&SenseInfo);
    SenseInfo += 10; /* "view_mode " */

    switch (SenseInfo[0])
    {
    case 'h':
        Mem->ViewQuality = VQ_High;
        break; /* high */
    case 'l':
        Mem->ViewQuality = VQ_Low;
        break; /* low  */
    default:
        my_error("Unknown view quality");
    }

    Mem->LastViewWidth = Mem->ViewWidth;
    Mem->ViewWidthTime = time;
    get_next_word(&SenseInfo);
    switch (SenseInfo[1])
    {
    case 'o':
        Mem->ViewWidth = VW_Normal;
        break; /* normal */
    case 'a':
        Mem->ViewWidth = VW_Narrow;
        break; /* narrow */
    case 'i':
        Mem->ViewWidth = VW_Wide;
        break; /* wide   */
    default:
        my_error("Unknown view quality");
    }

    float stamina = get_float(&SenseInfo);
    float effort = get_float(&SenseInfo);
    float speed = get_float(&SenseInfo);
    float head_angle = get_float(&SenseInfo);

    int kicks = get_int(&SenseInfo);
    int dashes = get_int(&SenseInfo);
    int turns = get_int(&SenseInfo);
    int says = get_int(&SenseInfo);
    int turn_necks = get_int(&SenseInfo);

    Mem->SetMySensedInfo(stamina, effort, speed, head_angle, kicks, dashes, turns, says, turn_necks, time);
}

/****************************************************************************************/

#define NOCHNGINFO -500
#define NOFACEINFO -500

void Parse_Sight(Time time, char *SightInfo)
{
    float dist, ang;
    float dirChng;
    float distChng;
    ObjType object_type;
    char player_side;
    Unum player_number;
    float facedir;
    float neckdir;
    MarkerType marker;
    SideLine line;
    Vqual view_qual;
    MarkerType closestMarker = No_Marker;
    Bool processThisMarker;
    float closestMarkerDist;
    /* float motionInfoDist = 1000; */

    while (*SightInfo != ')')
    {

        dirChng = NOCHNGINFO;
        facedir = NOFACEINFO;
        neckdir = NOFACEINFO;
        player_number = player_side = 0;

        get_word(&SightInfo); /* " ((" */

        if (*SightInfo == 'g')
        {
            object_type = OBJ_Marker;
            SightInfo += 5; /* "goal " */
            if (*SightInfo == 'r')
                marker = Goal_R;
            else if (*SightInfo == 'l')
                marker = Goal_L;
            else
                my_error("goal ?");
        }
        else if (*SightInfo == 'G')
        {
            object_type = OBJ_Marker_Behind;
            marker = Mem->ClosestGoal();
        }
        else if (*SightInfo == 'f')
        {
            object_type = OBJ_Marker;
            SightInfo += 5; /* "flag " */
            if (*SightInfo == 'r')
            {
                SightInfo += 2;
                if (*SightInfo == '0')
                    marker = Flag_R0;
                else if (*SightInfo == 'b')
                {
                    SightInfo += 1;
                    if (*SightInfo == ')')
                        marker = Flag_RB;
                    else
                    {
                        SightInfo += 1;
                        if (*SightInfo == '1')
                            marker = Flag_RB10;
                        else if (*SightInfo == '2')
                            marker = Flag_RB20;
                        else if (*SightInfo == '3')
                            marker = Flag_RB30;
                        else
                            my_error("flag r b ?");
                    }
                }
                else if (*SightInfo == 't')
                {
                    SightInfo += 1;
                    if (*SightInfo == ')')
                        marker = Flag_RT;
                    else
                    {
                        SightInfo += 1;
                        if (*SightInfo == '1')
                            marker = Flag_RT10;
                        else if (*SightInfo == '2')
                            marker = Flag_RT20;
                        else if (*SightInfo == '3')
                            marker = Flag_RT30;
                        else
                            my_error("flag r t ?");
                    }
                }
                else
                    my_error("flag r ?");
            }
            else if (*SightInfo == 'l')
            {
                SightInfo += 2;
                if (*SightInfo == '0')
                    marker = Flag_L0;
                else if (*SightInfo == 'b')
                {
                    SightInfo += 1;
                    if (*SightInfo == ')')
                        marker = Flag_LB;
                    else
                    {
                        SightInfo += 1;
                        if (*SightInfo == '1')
                            marker = Flag_LB10;
                        else if (*SightInfo == '2')
                            marker = Flag_LB20;
                        else if (*SightInfo == '3')
                            marker = Flag_LB30;
                        else
                            my_error("flag l b ?");
                    }
                }
                else if (*SightInfo == 't')
                {
                    SightInfo += 1;
                    if (*SightInfo == ')')
                        marker = Flag_LT;
                    else
                    {
                        SightInfo += 1;
                        if (*SightInfo == '1')
                            marker = Flag_LT10;
                        else if (*SightInfo == '2')
                            marker = Flag_LT20;
                        else if (*SightInfo == '3')
                            marker = Flag_LT30;
                        else
                            my_error("flag l t ?");
                    }
                }
                else
                    my_error("flag l ?");
            }
            else if (*SightInfo == 't')
            {
                SightInfo += 2;
                if (*SightInfo == '0')
                    marker = Flag_T0;
                else if (*SightInfo == 'l')
                {
                    SightInfo += 2;
                    if (*SightInfo == '1')
                        marker = Flag_TL10;
                    else if (*SightInfo == '2')
                        marker = Flag_TL20;
                    else if (*SightInfo == '3')
                        marker = Flag_TL30;
                    else if (*SightInfo == '4')
                        marker = Flag_TL40;
                    else if (*SightInfo == '5')
                        marker = Flag_TL50;
                    else
                        my_error("flag t l ?");
                }
                else if (*SightInfo == 'r')
                {
                    SightInfo += 2;
                    if (*SightInfo == '1')
                        marker = Flag_TR10;
                    else if (*SightInfo == '2')
                        marker = Flag_TR20;
                    else if (*SightInfo == '3')
                        marker = Flag_TR30;
                    else if (*SightInfo == '4')
                        marker = Flag_TR40;
                    else if (*SightInfo == '5')
                        marker = Flag_TR50;
                    else
                        my_error("flag t r ?");
                }
                else
                    my_error("flag t ?");
            }
            else if (*SightInfo == 'b')
            {
                SightInfo += 2;
                if (*SightInfo == '0')
                    marker = Flag_B0;
                else if (*SightInfo == 'l')
                {
                    SightInfo += 2;
                    if (*SightInfo == '1')
                        marker = Flag_BL10;
                    else if (*SightInfo == '2')
                        marker = Flag_BL20;
                    else if (*SightInfo == '3')
                        marker = Flag_BL30;
                    else if (*SightInfo == '4')
                        marker = Flag_BL40;
                    else if (*SightInfo == '5')
                        marker = Flag_BL50;
                    else
                        my_error("flag b l ?");
                }
                else if (*SightInfo == 'r')
                {
                    SightInfo += 2;
                    if (*SightInfo == '1')
                        marker = Flag_BR10;
                    else if (*SightInfo == '2')
                        marker = Flag_BR20;
                    else if (*SightInfo == '3')
                        marker = Flag_BR30;
                    else if (*SightInfo == '4')
                        marker = Flag_BR40;
                    else if (*SightInfo == '5')
                        marker = Flag_BR50;
                    else
                        my_error("flag b r ?");
                }
                else
                    my_error("flag b ?");
            }
            else if (*SightInfo == 'c')
            {
                SightInfo += 1;
                if (*SightInfo == ')')
                    marker = Flag_C;
                else
                {
                    SightInfo += 1;
                    if (*SightInfo == 'b')
                        marker = Flag_CB;
                    else if (*SightInfo == 't')
                        marker = Flag_CT;
                    else
                        my_error("flag c ?");
                }
            }
            else if (*SightInfo == 'p')
            {
                SightInfo += 2;
                if (*SightInfo == 'r')
                {
                    SightInfo += 2;
                    if (*SightInfo == 't')
                        marker = Flag_PRT;
                    else if (*SightInfo == 'c')
                        marker = Flag_PRC;
                    else if (*SightInfo == 'b')
                        marker = Flag_PRB;
                    else
                        my_error("flag p r ?");
                }
                else if (*SightInfo == 'l')
                {
                    SightInfo += 2;
                    if (*SightInfo == 't')
                        marker = Flag_PLT;
                    else if (*SightInfo == 'c')
                        marker = Flag_PLC;
                    else if (*SightInfo == 'b')
                        marker = Flag_PLB;
                    else
                        my_error("flag p l ?");
                }
                else
                    my_error("flag p ?");
            }
            else if (*SightInfo == 'g')
            {
                SightInfo += 2;
                if (*SightInfo == 'l')
                {
                    SightInfo += 2;
                    if (*SightInfo == 't')
                        marker = Flag_GLT;
                    else if (*SightInfo == 'b')
                        marker = Flag_GLB;
                    else
                        my_error("flag g l ?");
                }
                else if (*SightInfo == 'r')
                {
                    SightInfo += 2;
                    if (*SightInfo == 't')
                        marker = Flag_GRT;
                    else if (*SightInfo == 'b')
                        marker = Flag_GRB;
                    else
                        my_error("flag g r ?");
                }
                else
                    my_error("flag g ?");
            }
            else
                my_error("flag ?");
        }
        else if (*SightInfo == 'F')
        {
            object_type = OBJ_Marker_Behind;
            marker = Mem->ClosestFlagTo(); /* could be No_Marker */
        }
        else if (*SightInfo == 'l')
        {
            object_type = OBJ_Line;
            SightInfo += 5; /* "line " */
            if (*SightInfo == 'r')
                line = SL_Right;
            else if (*SightInfo == 'l')
                line = SL_Left;
            else if (*SightInfo == 't')
                line = SL_Top;
            else if (*SightInfo == 'b')
                line = SL_Bottom;
            else
                my_error("line ?");
        }
        else if (*SightInfo == 'p' || *SightInfo == 'P')
        {
            object_type = OBJ_Player;
            SightInfo += 6; /* "player" */
            if (*SightInfo == ' ')
            { /* there's a team */
                SightInfo++;
                if (!strncmp(SightInfo, Mem->MyTeamName, Mem->MyTeamNameLen))
                    player_side = Mem->MySide;
                else
                {
                    if (Mem->TheirTeamName[0] == '\n')
                    {
                        int a = 0;
                        while (isalpha(*SightInfo))
                            Mem->TheirTeamName[a++] = *SightInfo++;
                    }
                    player_side = Mem->TheirSide;
                }
                while (*SightInfo != ' ' && *SightInfo != ')')
                    SightInfo++; /* advance past team name */
                if (*SightInfo == ' ')
                { /* there's a number */
                    player_number = get_int(&SightInfo);
                }
            }
        }
        else if (*SightInfo == 'b' || *SightInfo == 'B')
            object_type = OBJ_Ball;
        else
            my_error("unknown object");

        advance_to(')', &SightInfo); /* advance to end of object */

        /************************************/

        ang = get_float(&SightInfo);

        if (*SightInfo != ')')
        { /* 'high' quality     */
            view_qual = VQ_High;
            dist = ang;
            ang = get_float(&SightInfo);
        }
        else
        {
            printf("%s", SightInfo - 30);
            view_qual = VQ_Low;
        }

        if (view_qual != Mem->ViewQuality)
            my_error("View quality %d correct?", view_qual);

        if (*SightInfo != ')')
        {
            distChng = get_float(&SightInfo);
            dirChng = get_float(&SightInfo);
        }

        if (*SightInfo != ')')
        {
            if (object_type != OBJ_Player)
                my_error("Only players should have facedir");
            facedir = get_float(&SightInfo);
            neckdir = get_float(&SightInfo);
        }

        if (*SightInfo != ')')
            my_error("Should be done with object info here");
        SightInfo++; /* ")" */

        /************************************/

        switch (object_type)
        {
        case OBJ_Marker:
        case OBJ_Marker_Behind:
            /* Want to save 2 closest for triangulation  */
            /* don't want marker_behind unless necessary */

            /* If it was a Marker_Behind and we don't know which one */
            if (marker == No_Marker)
            {
                if (object_type != OBJ_Marker_Behind)
                    my_error("Should know the marker");
                break;
            }

            processThisMarker = FALSE;
            if (view_qual == VQ_Low)
            {                                  /* Low quality   */
                                               /* DON'T BOTHER PROCESSING ANY??? I don't think it helps ... */
                                               /* COULD process 2---then triangulate */
                                               /*if ( closestMarkerDist > 0 ){ */
                /*  closestMarkerDist = 0;  */ /* Only process 1*/
                                               /*  processThisMarker = TRUE; */
                                               /*}*/
            }
            else
            { /* high quality  */
                if (closestMarker == No_Marker || dist < closestMarkerDist)
                {
                    closestMarker = marker;
                    closestMarkerDist = dist;
                    processThisMarker = TRUE;
                    Mem->ClosestMarker = marker;
                }
                /* Don't bother with marker motion info -- get it from sense_body and my angle
                if ( dirChng != NOCHNGINFO && dist < motionInfoDist ){
                  motionInfoDist = dist;
                  processThisMarker = TRUE;
                  Mem->ClosestMotionMarker = marker;
                }
                */
            }
            if (processThisMarker)
            {
                if (view_qual == VQ_Low) /* low quality   */
                    Mem->SeeMarker(marker, ang, time);
                else /* if (dirChng == NOCHNGINFO) */        /* high quality  */
                    Mem->SeeMarker(marker, dist, ang, time); /* No motion info*/
                                                             /* else
                                                               Mem->SeeMarker(marker, dist, ang, distChng, dirChng, time); */
            }
            break;
        case OBJ_Line:
            if (*SightInfo != ')')
                /* There's another line coming.  Assuming lines happen
                   last in the visual string and the closer line comes first */
                ;
            else if (view_qual == VQ_Low) /* low quality   */
                Mem->SeeLine(line, ang, time);
            else /* high quality  */
                Mem->SeeLine(line, dist, ang, time);
            break;
        case OBJ_Ball:
            if (view_qual == VQ_Low) /* low quality   */
                Mem->SeeBall(ang, time);
            else if (dirChng == NOCHNGINFO) /* high quality  */
                Mem->SeeBall(dist, ang, time);
            else /* know direction*/
                Mem->SeeBall(dist, ang, distChng, dirChng, time);
            break;
        case OBJ_Player:
            if (!player_side)
            {                            /* Too far for team or num */
                if (view_qual == VQ_Low) /* low quality   */
                    Mem->SeePlayer(ang, time);
                else if (dirChng == NOCHNGINFO) /* high quality  */
                    Mem->SeePlayer(dist, ang, time);
                else /* know direction*/
                    my_error("Shouldn't know dirChng when the player's far");
            }

            else
            {
                if (!player_number)
                {                            /* Too far for number     */
                    if (view_qual == VQ_Low) /* low quality   */
                        Mem->SeePlayer(player_side, ang, time);
                    else if (dirChng == NOCHNGINFO) /* high quality  */
                        Mem->SeePlayer(player_side, dist, ang, time);
                    else /* know direction*/
                        my_error("Shouldn't know dirChng when the team member's far");
                }

                else
                {                            /* Know side AND number   */
                    if (view_qual == VQ_Low) /* low quality   */
                        Mem->SeePlayer(player_side, player_number, ang, time);
                    else if (dirChng == NOCHNGINFO)
                    { /* high quality  */
                        printf("%s\n", SightInfo - 30);
                        my_error("Should know dirChng when know number");
                        Mem->SeePlayer(player_side, player_number, dist, ang, time);
                    }
                    else /* know direction*/
                        Mem->SeePlayer(player_side, player_number, dist, ang, distChng, dirChng, facedir, neckdir, time);
                }
            }
        }
    }
}

/****************************************************************************************/

void Parse_Sound(Time time, char *SoundInfo)
{
    if (SoundInfo[1] == 'r')
    {                   /* Referee or Coach message */
        SoundInfo += 9; /* " referee " */
        if (strncmp(SoundInfo, "training", 8) == 0)
            Parse_Trainer_Sound(SoundInfo);
        else if (islower(SoundInfo[0]))
            Parse_Referee_Sound(SoundInfo);
        else
            my_error("Referee sounds should start with lower case letters!");
        return;
    }
    else if (SoundInfo[1] == 'o')
    {                    /* Online coach message */
        SoundInfo += 14; /* online_coach_ */
        if (SoundInfo[0] == Mem->MySide)
        {
            advance_to(' ', &SoundInfo); /* advance to end of side */
            SoundInfo++;
            Parse_My_Coach_Sound(time, SoundInfo);
        }
        else if (SoundInfo[0] == Mem->TheirSide)
        {
            advance_to(' ', &SoundInfo); /* advance to end of side */
            SoundInfo++;
            Parse_Their_Coach_Sound(time, SoundInfo);
        }
        else
            my_error("online_coach_?");
        return;
    }

    /* Here you could do something with the heard message */
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void Parse_Referee_Sound(char *msg)
{
    std::cout << msg << std::endl;
    // printf("Parsing ref sound: '%s\'\n", msg);
    msg[strlen(msg) - 1] = 0; /* cut off the newline and paren */
    Mem->LogAction3(200, "Referee message: %s", msg);
    switch (msg[0])
    {
    case 'p':
        Mem->SetPlayMode(PM_Play_On);
        break; /* play_on */
    case 'k':
        if (msg[5] == 'i')
        { /* kick_in */
            if (msg[8] == Mem->MySide)
                Mem->SetPlayMode(PM_My_Kick_In);
            else if (msg[8] == Mem->TheirSide)
                Mem->SetPlayMode(PM_Their_Kick_In);
            else
                my_error("kick_in_?");
        }
        else if (msg[5] == 'o')
        { /* kick_off */
            if (msg[9] == Mem->MySide)
                Mem->SetPlayMode(PM_My_Kick_Off);
            else if (msg[9] == Mem->TheirSide)
                Mem->SetPlayMode(PM_Their_Kick_Off);
            else
                my_error("kick_off_?");
        }
        else
            my_error("referee k..?");
        break;
    case 'g':
        if (msg[5] == 'k')
        { /* goal_kick */
            if (msg[10] == Mem->MySide)
                Mem->SetPlayMode(PM_My_Goal_Kick);
            else if (msg[10] == Mem->TheirSide)
                Mem->SetPlayMode(PM_Their_Goal_Kick);
            else
                my_error("goal_kick_?");
        }
        else if (msg[5] == 'e')
        { /* goalie_catch_ball */
            if (msg[18] == Mem->MySide)
                Mem->SetPlayMode(PM_My_Goalie_Free_Kick);
            else if (msg[18] == Mem->TheirSide)
                Mem->SetPlayMode(PM_Their_Goalie_Free_Kick);
            else
                my_error("goalie_catch_ball_?");
        }
        else if (msg[5] == Mem->MySide)
        { /* goal */
            Mem->MyScore++;
            // Mem->MyScore = get_int(&msg[7]);
            to_ori_pos();
            Mem->KickOffMode = KO_Theirs;
            Mem->SetPlayMode(PM_Before_Kick_Off);
        }
        else if (msg[5] == Mem->TheirSide)
        {
            Mem->TheirScore++;
            // Mem->TheirScore = get_int(&msg[7]);
            to_ori_pos();
            Mem->KickOffMode = KO_Mine;
            Mem->SetPlayMode(PM_Before_Kick_Off);
        }
        else
            my_error("referee g..?");
        break;
    case 'c': /* corner_kick */
        if (msg[12] == Mem->MySide)
            Mem->SetPlayMode(PM_My_Corner_Kick);
        else if (msg[12] == Mem->TheirSide)
            Mem->SetPlayMode(PM_Their_Corner_Kick);
        else
            my_error("corner_kick_?");
        break;
    case 'd':
        Mem->SetPlayMode(PM_Drop_Ball);
        break; /* drop_ball */
    case 'o':  /* offside */
        if (msg[8] == Mem->MySide)
            Mem->SetPlayMode(PM_Their_Offside_Kick);
        else if (msg[8] == Mem->TheirSide)
            Mem->SetPlayMode(PM_My_Offside_Kick);
        else
            my_error("offside_?");
        break;
    case 'f':
        if (msg[5] == 'k')
        { /* free_kick */
            if (msg[10] == Mem->MySide)
                Mem->SetPlayMode(PM_My_Free_Kick);
            else if (msg[10] == Mem->TheirSide)
                Mem->SetPlayMode(PM_Their_Free_Kick);
            else
                my_error("free_kick_?");
        }
        else if (msg[5] == Mem->MySide) /* foul */
            ;
        else if (msg[5] == Mem->TheirSide)
            ;
        else
            my_error("referee f..?");
        break;
    case 'h':                           /* half_time */
        Mem->SetPlayMode(PM_Half_Time); /* play_mode to before_kick_off        */
        if (Mem->MySide == 'l')
            Mem->KickOffMode = KO_Theirs;
        else
            Mem->KickOffMode = KO_Mine;
        break;
    case 'b':
        Mem->SetPlayMode(PM_Before_Kick_Off);
        break; /* before_kick_off */
    case 't':
        if (msg[5] == 'u')
        { /* time_up */
            Mem->SetPlayMode(PM_Time_Up);
        }
        else if (msg[5] == 'o') /* time_over */
        {
            break;
        }
        else if (msg[5] == 'e')
        { /* time_extended */
            Mem->SetPlayMode(PM_Extended_Time);
            if (Mem->MySide == 'l')
                Mem->KickOffMode = KO_Mine;
            else
                Mem->KickOffMode = KO_Theirs;
        }
        else
            my_error("referee t..?");
        break;
    default:
        my_error("Referee msg ????");
    }
}

/****************************************************************************************/

/* the trainer send a string that is essentially command line options */
void Parse_Trainer_Sound(char *msg)
{
    msg += 9;                 /* 'training ' */
    msg[strlen(msg) - 1] = 0; /* cut off the newline and paren */
    Mem->GetOptions(msg);
    Mem->GetBall()->forget(); /* set 0 confidence in ball pos */
    std::cout << "Incorp trainer message" << std::endl;
    Mem->LogAction2(175, "Incorporated trainer sound");
}

/****************************************************************************************/

void Parse_My_Coach_Sound(Time time, char *msg)
{
    /* Here you can parse the coach message */
}

/****************************************************************************************/

void Parse_Their_Coach_Sound(Time time, char *msg)
{
    /* Here you can parse their coach's messages */
}

/* -*- Mode: C++ -*- */

/* test.C
 * CMUnited99 (soccer client for Robocup99)
 * Peter Stone <pstone@cs.cmu.edu>
 * Computer Science Department
 * Carnegie Mellon University
 * Copyright (C) 1999 Peter Stone
 *
 * CMUnited-99 was created by Peter Stone, Patrick Riley, and Manuela Veloso
 *
 * You may copy and distribute this program freely as long as you retain this notice.
 * If you make any changes or have any comments we would appreciate a message.
 * For more information, please see http://www.cs.cmu.edu/~robosoccer/
 */

/*****************************************************************************************/

void test_scan_with_body()
{
    if (Mem->PlayMode == PM_Before_Kick_Off)
    {
        move(-1, 0);
        my_stamp;
        printf("%.1f\n", Mem->MyBodyAng());
        return;
    }

    if (Mem->TimeToTurnForScan() && Mem->CurrentTime < 50)
    {
        my_stamp;
        printf("turning.  Last sight at %d\n", Mem->LastSightTime.t);
        turn(Mem->MyViewAngle() * 2);
    }
    my_stamp;
    printf("        %.1f\n", Mem->MyBodyAng());
}

/*****************************************************************************************/

void test_random_movement_in_rectangle(Rectangle *rect)
{
    static int First_time = TRUE;
    if (First_time)
    {
        Vector point = rect->random();
        point = Mem->PositionToKickoffPosition(point);
        move(point.x, point.y);
        First_time = FALSE;
        return;
    }

    if (!Mem->MyConf())
    {
        turn(30);
        return;
    }

    if (rect->DistanceToEdge(Mem->MyPos()) < 5)
        test_go_to_point(rect->Center(), 0);
    else if (Mem->MyStamina() >= Mem->EffortDecThreshold)
        test_random_movement();
}

/*****************************************************************************************/

void test_1v1()
{
    static int First_time = TRUE;
    if (First_time)
    {
        if (Mem->PlayMode == PM_Before_Kick_Off)
            move(-10, 0);
        First_time = FALSE;
        return;
    }

    if (Mem->BallPositionValid())
    {
        if (Mem->BallKickable())
        {
            hold_ball();
        }
        else
        {
            get_ball();
        }
    }
    else
        face_neck_and_body_to_ball();
}

/*****************************************************************************************/

void test_volley()
{
    static int First_time = TRUE;
    if (First_time)
    {
        move(-20, 0);
        First_time = FALSE;
        return;
    }

    if (Mem->BallPositionValid() && Mem->BallX() < -1)
        test_go_to_ball();
    else
        test_go_to_point(Vector(-20, 0), 3);
}

/*****************************************************************************************/

void test_go_to_ball(AngleDeg kick_angle)
{
    if (Mem->BallKickable())
    {
        if (!Mem->TestVersion)
        {
            if (!smart_kick_hard(kick_angle, KM_QuickestRelease))
                my_error("test_go_to_ball: kick failed\n");
        }
        else
            kick(100, kick_angle);
    }
    else if (Mem->MyConf() && Mem->BallPositionValid())
    {
        if (Mem->TestVersion)
            test_go_to_point(Mem->BallAbsolutePosition(), Mem->SP_kickable_area);
        else
        {
            if (Mem->MyInterceptionAble())
                test_go_to_point(Mem->MyInterceptionPoint(), 0);
            else
                test_face_ball();
        }
    }
    else
        turn(90);
}

/*****************************************************************************************/

void test_go_to_ball()
{
    if (Mem->MyConf())
        test_go_to_ball(Mem->MarkerAngleFromBody(Mem->RM_Their_Goal));
    else
        turn(90);
}

/*****************************************************************************************/

void test_go_to_point(Vector p, float buffer, float dash_power)
{
    if (!Mem->MyConf())
    {
        turn(Mem->MyViewAngle());
        return;
    }

    if (Mem->DistanceTo(p) < buffer)
    {
        // turn(30);
        test_face_ball();
        return;
    }

    float target_ang = Mem->AngleToFromBody(p);
    float target_dist = Mem->DistanceTo(p);

    if (1 && !Mem->ClockStopped)
    { /* dodge players */
        PlayerObject *player;
        float dodge_dist = Min(Mem->CP_dodge_distance_buffer, target_dist);
        AngleDeg dodge_ang = Mem->CP_dodge_angle_buffer;
        if ((player = Mem->GetPlayerWithin(dodge_dist, dodge_ang, 0, target_ang - dodge_ang)) != NULL)
        {
            if (Mem->NumPlayersWithin(dodge_dist, 2 * dodge_ang))
            {
                /* Target close, so no players will be within in the next iteration ==> dash */
                Vector new_target = Mem->BodyPolar2Gpos(Mem->SP_player_size, player->get_ang_from_body() + 90);
                test_go_to_point(new_target, 0);
            }
            else
            {
                dash(Mem->CP_dodge_power);
            }
            return;
        }
        if ((player = Mem->GetPlayerWithin(dodge_dist, dodge_ang, 0, target_ang + dodge_ang)) != NULL)
        {
            if (Mem->NumPlayersWithin(dodge_dist, 2 * dodge_ang))
            {
                /* Target close, so no players will be within in the next iteration ==> dash */
                Vector new_target = Mem->BodyPolar2Gpos(Mem->SP_player_size, player->get_ang_from_body() - 90);
                test_go_to_point(new_target, 0);
            }
            else
            {
                dash(Mem->CP_dodge_power);
            }
            return;
        }
    }

    if (fabs(target_ang) > Mem->CP_max_go_to_point_angle_err)
    {
        turn(target_ang);
        return;
    }

    if (Mem->MyStamina() >= Mem->EffortDecThreshold)
        dash(Min(dash_power, Mem->MyStamina() - Mem->EffortDecThreshold));
    else
    {
        my_stamp;
        printf("recovering\n");
    } // turn(180);
}

/*****************************************************************************************/

void test_face_ball()
{
    if (Mem->BallPositionValid())
        turn(Mem->BallAngleFromNeck());
    else
        turn(Mem->MyViewAngle());
}

/*****************************************************************************************/

void test_random_movement()
{
    static int First_time = FALSE;
    if (First_time)
    {
        move(range_random(-50, 0), range_random(-Mem->SP_pitch_width / 2, Mem->SP_pitch_width / 2));
        First_time = FALSE;
        return;
    }

    /* if      ( !int_random(100) ) change_view(VW_Wide);
       else if ( !int_random(100) ) change_view(VW_Normal);
       else if ( !int_random(100) ) change_view(VW_Narrow); */

    if (Mem->ClockStopped)
        return;

    if (Mem->BallKickable())
    {
        kick(range_random(0, 100), Mem->BallAngleFromBody());
        return;
    }

    if (Mem->MyConf() && Mem->BallPositionValid() && Mem->BallDistance() < 9)
    {
        if (fabs(Mem->BallAngleFromBody()) > 1)
            turn(Mem->BallAngleFromBody());
        else
            dash(100);
        return;
    }

    if (int_random(4))
        dash(range_random(20, 30));
    else
        turn(range_random(-90, 90));

    return;
}

/*****************************************************************************************/

void test_straight_to_ball()
{
    static int First_time = TRUE;
    if (First_time)
    {
        if (Mem->MySide == 'l')
            turn(90);
        else
            turn(-90);
        First_time = FALSE;
        return;
    }

    if (Mem->CurrentTime.t && Mem->CurrentTime.t < 10 && Mem->GetBall()->pos_valid())
    {
        turn(Mem->GetBall()->get_ang_from_body());
    }
    else
        dash(Mem->CurrentTime.t % 100);
    return;
}

/*****************************************************************************************/

void test_run_straight()
{
    static Bool GO = TRUE;

    if (!(Mem->CurrentTime.t % 200))
        GO = TRUE;

    if (Mem->MyStamina() > 600 && GO)
        dash(60);
    else if (Mem->CurrentTime.t)
    {
        GO = FALSE;
        printf("stamina = %f, effort = %f\n", Mem->MyStamina(), Mem->MyEffort());
    }
}

/*****************************************************************************************/

void test_turn_and_dash_slow()
{
    static int First_time = TRUE;
    if (First_time)
    {
        if (Mem->MySide == 'l')
            turn(90);
        else
            turn(-90);
        First_time = FALSE;
        return;
    }

    if (Mem->CurrentTime.t % 5 || !Mem->CurrentTime.t)
        return;

    if (!(Mem->CurrentTime.t % 50))
        turn(40);
    else
        dash(100);

    return;
}

/*****************************************************************************************/

void test_print_ball()
{
    my_stamp;

    if (Mem->BallPositionValid())
    {
        printf(" Ball at %.1f %.1f (%.2f)", Mem->BallX(), Mem->BallY(), Mem->BallPositionValid());
        if (Mem->BallVelocityValid())
            printf("     Speed %.2f dir %.1f (%.2f)",
                   Mem->BallSpeed(), Mem->BallAbsoluteHeading(), Mem->BallVelocityValid());
        else
            printf("    velocity not valid");
    }
    else
        printf(" Ball not valid");

    printf("\n");
}

/*****************************************************************************************/

void test_print_positions()
{
    Bool me = TRUE;
    Bool ball = FALSE;
    Bool players = TRUE;

    if (me)
    {
        my_stamp;
        printf("at %d I'm at (%.1f %.1f) facing %.1f with velocity (%f %.1f)\n",
               Mem->CurrentTime.t, Mem->MyX(), Mem->MyY(), Mem->MyBodyAng(), Mem->MySpeed(), Mem->MyDir());
    }

    if (ball)
    {
        if (Mem->GetBall()->pos_valid())
            printf("Ball: (%f %f)       ", Mem->GetBall()->get_x(), Mem->GetBall()->get_y());

        if (Mem->GetBall()->vel_valid())
            printf("velocity: (%f %f)\n", Mem->GetBall()->get_speed(), Mem->GetBall()->get_abs_heading());
        else if (Mem->GetBall()->pos_valid())
            printf("\n");
    }

    if (players)
    {
        for (int i = 0; i < Mem->NumTeammates(); i++)
        {
            my_stamp;
            printf("T %d at %.1f %.1f facing %.1f\n", Mem->Teammates()[i]->unum, Mem->Teammates()[i]->get_x(), Mem->Teammates()[i]->get_y(), Mem->Teammates()[i]->get_abs_body_ang());
        }
        for (int i = 0; i < Mem->NumOpponents(); i++)
        {
            my_stamp;
            printf("O %d at %.1f %.1f facing %.1f\n", Mem->Opponents()[i]->unum, Mem->Opponents()[i]->get_x(), Mem->Opponents()[i]->get_y(), Mem->Opponents()[i]->get_abs_body_ang());
        }
        for (int i = 0; i < Mem->NumTeamlessPlayers(); i++)
        {
            my_stamp;
            printf("? %d at %.1f %.1f facing %.1f\n", Mem->TeamlessPlayers()[i]->unum, Mem->TeamlessPlayers()[i]->get_x(), Mem->TeamlessPlayers()[i]->get_y(), Mem->TeamlessPlayers()[i]->get_abs_body_ang());
        }
    }
}

/*****************************************************************************************/

void test_time()
{
    Time tmp = Mem->CurrentTime - 2;
    printf("time: %d %d         %d %d      (%d %d)",
           Mem->CurrentTime.t, Mem->CurrentTime.s, tmp.t, tmp.s, Mem->LastStartClockTime.t, Mem->LastStartClockTime.s);
    if (Mem->CurrentTime.t > 5)
    {
        int tmp2 = Mem->CurrentTime - Time(5, 0);
        printf(" --------  %d)", tmp2);
    }
    printf("\n");
}

/***************************************************************************/
/* Pat added for parameter tuning */
/**************************************************************************/
void test_turnball()
{
    if (!Mem->MyConf() || !Mem->BallPositionValid())
        scan_field_with_body();
    else if (Mem->BallKickable())
    {
        TurnballTo(Mem->BallAngleFromBody() + 180, TURN_CCW);
    }
    else
        ; /* if we can't kick it, do nothing */
}

void test_turnball2()
{
    static TurnDir turn_dir = TURN_CW;
    if (!Mem->MyConf() || !Mem->BallPositionValid())
    {
        /* switch the turn direction when we lose the ball,
           this effectively randomizes it */
        turn_dir = (turn_dir == TURN_CW) ? TURN_CCW : TURN_CW;
        scan_field_with_body();
    }
    else if (Mem->BallKickable())
    {
        TurnballTo(-Mem->MyBodyAng(), turn_dir);
    }
    else
        ; /* if we can't kick it, do nothing */
}

#ifdef DEBUG_OUTPUT
#define DEBUG_TEST_KICK(x)
#else
#define DEBUG_TEST_KICK(x)
#endif
void test_hard_kick(KickMode km)
{
    static float bally;
    static Bool LastActionKick = FALSE;

    if (!Mem->MyConf())
    {
        DEBUG_TEST_KICK(cout << "Time: " << Mem->CurrentTime.t << "\tScanning field" << endl);
        scan_field_with_body();
        LastActionKick = FALSE;
    }
    else if (Mem->PlayMode == PM_Before_Kick_Off)
    {
        move(-10, 0);
    }
    else
    {
        DEBUG_TEST_KICK(cout << "Time: " << Mem->CurrentTime.t << "\tMyAng: " << Mem->MyAng() << endl);
        if (fabs(GetNormalizeAngleDeg(Mem->MyBodyAng() -
                                      signf(bally) *
                                          Mem->CP_hardest_kick_player_ang)) > 9 &&
            (!Mem->BallKickable() || Mem->BallVelocityValid()) &&
            !LastActionKick)
        {
            DEBUG_TEST_KICK(cout << " Turning: target: "
                                 << signf(bally) * Mem->CP_hardest_kick_player_ang
                                 << "\tturning: "
                                 << signf(bally) * Mem->CP_hardest_kick_player_ang - Mem->MyAng() << endl);
            turn(signf(bally) * Mem->CP_hardest_kick_player_ang - Mem->MyBodyAng());
            LastActionKick = FALSE;
        }
        else if (Mem->BallKickable())
        {
            DEBUG_TEST_KICK(cout << "calling smart_kick_hard" << endl);
            smart_kick_hard_abs(0.0, km);
            if (Mem->Action->type == CMD_kick)
                LastActionKick = TRUE;
            else
                LastActionKick = FALSE;
            bally = Mem->BallY();
        }
        else
            LastActionKick = FALSE; /* if we can't kick or turn, do nothing */
    }
}

void test_intercept()
{
    if (!Mem->MyConf())
    {
        Mem->LogAction2(10, "Lost myself, scanning field");
        scan_field_with_body();
    }
    else if (!Mem->BallPositionValid())
    {
        Mem->LogAction2(10, "Lost the ball, scanning field");
        scan_field_with_body();
    }
    else if (Mem->BallKickable())
    {
        Mem->LogAction2(10, "Have ball, holding onto it");
        hold_ball();
    }
    else
    {
        Mem->LogAction2(10, "Chasing the ball");
        get_ball();
    }
}

void test_go_to_static_ball()
{
    if (!Mem->MyConf() || !Mem->BallPositionValid())
        scan_field_with_body();
    else if (Mem->BallKickable() && Mem->LastActionType() == CMD_kick)
    {
        smart_kick_hard_abs((-Mem->BallAbsolutePosition()).dir(), KM_HardestKick);
    }
    else
    {
        if (go_to_static_ball((-Mem->BallAbsolutePosition()).dir()))
            smart_kick_hard_abs((-Mem->BallAbsolutePosition()).dir(), KM_HardestKick);
        // cout << Mem->CurrentTime.t << " I'm ready to kick" << endl;
    }
}

void test_pred_cycles_to_point()
{
    static Vector pt(0, 0);
    static Rectangle r(Vector(0, 0), Vector(20, 20));

    if (!Mem->MyConf())
        scan_field_with_body();
    else
    {
        std::cout << "Time: " << std::setw(4) << Mem->CurrentTime.t << '\t'
                  << "predict " << std::setw(3) << Mem->PredictedCyclesToPoint(pt) << " cycles "
                  << "to point " << pt << std::endl;
        if (go_to_point(pt, Mem->CP_at_point_buffer) == AQ_ActionNotQueued)
        {
            std::cout << "I got there!" << std::endl;
            pt = r.random();
        }
    }
}

void test_log_action()
{
    if (!Mem->MyConf())
    {
        Mem->LogAction2(10, "I'm lost! Looking around");
        scan_field_with_body();
    }
    else if (!Mem->BallPositionValid())
    {
        Mem->LogAction2(10, "I lost the ball! Looking around");
        scan_field_with_body();
    }
    else if (Mem->BallKickable())
    {
        float kick_angle;
        Mem->LogAction2(10, "Ball is kickable!");
        if (Mem->MyX() > 0)
        {
            Mem->LogAction2(30, "On offensive half, kicking to goal");
            kick_angle = Mem->MarkerAngleFromBody(Mem->RM_Their_Goal);
        }
        else
        {
            Mem->LogAction2(30, "On defensive half, kicking to corner");
            kick_angle = Mem->MarkerAngleFromBody(Mem->RM_RF_Flag);
        }

        if (!smart_kick_hard(kick_angle, KM_QuickestRelease))
            my_error("test_go_to_ball: kick failed\n");
    }
    else
    {
        Mem->LogAction2(10, "I see the ball! going for it");
        test_go_to_point(Mem->BallAbsolutePosition(), Mem->SP_kickable_area);
    }
}