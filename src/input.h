#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

#include "typedefs.h"


enum KeyboardKeyState {
    KeyboardKeyState_Depressed = 0, // No recent event, key is still up
    KeyboardKeyState_Released = 1,  // Last event was a released event
    KeyboardKeyState_Held = 2,      // No recent event, key is still down
    KeyboardKeyState_Pressed = 3,    // Last event was a pressed event
    KeyboardKeyState_COUNT = 4
};

struct KeyBoardState {
    enum KeyboardKeyState up;
    enum KeyboardKeyState down;
    enum KeyboardKeyState right;
    enum KeyboardKeyState left;
    enum KeyboardKeyState a;
    enum KeyboardKeyState s;
    enum KeyboardKeyState d;
    enum KeyboardKeyState w;
    enum KeyboardKeyState q;
    enum KeyboardKeyState e;
    enum KeyboardKeyState r;

    enum KeyboardKeyState one;
    enum KeyboardKeyState two;
    enum KeyboardKeyState three;
    enum KeyboardKeyState four;
    enum KeyboardKeyState five;
    enum KeyboardKeyState six;
    enum KeyboardKeyState seven;
    enum KeyboardKeyState eight;
};

void ClearKeyboardState(struct KeyBoardState* kbs);

void DecayKeyboardState(struct KeyBoardState* kbs);

bool IsPressed(enum KeyboardKeyState state);

bool IsNewlyPressed(enum KeyboardKeyState state);

#endif