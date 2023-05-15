#include "input.h"

#include <stddef.h>

void ClearKeyboardState(struct KeyBoardState* kbs) {
    kbs->up = KeyboardKeyState_Depressed;
    kbs->down = KeyboardKeyState_Depressed;
    kbs->right = KeyboardKeyState_Depressed;
    kbs->left = KeyboardKeyState_Depressed;
    kbs->a = KeyboardKeyState_Depressed;
    kbs->s = KeyboardKeyState_Depressed;
    kbs->d = KeyboardKeyState_Depressed;
    kbs->w = KeyboardKeyState_Depressed;
    kbs->q = KeyboardKeyState_Depressed;
    kbs->e = KeyboardKeyState_Depressed;
    kbs->r = KeyboardKeyState_Depressed;

    kbs->one = KeyboardKeyState_Depressed;
    kbs->two = KeyboardKeyState_Depressed;
    kbs->three = KeyboardKeyState_Depressed;
    kbs->four = KeyboardKeyState_Depressed;
    kbs->five = KeyboardKeyState_Depressed;
    kbs->six = KeyboardKeyState_Depressed;
    kbs->seven = KeyboardKeyState_Depressed;
    kbs->eight = KeyboardKeyState_Depressed;
}

void DecayKeyboardState(struct KeyBoardState* kbs) {
    static enum KeyboardKeyState to_depressed_state[KeyboardKeyState_COUNT] = {
        KeyboardKeyState_Depressed,
        KeyboardKeyState_Depressed,
        KeyboardKeyState_Held,
        KeyboardKeyState_Held
    };

    kbs->up = to_depressed_state[kbs->up];
    kbs->down =  to_depressed_state[kbs->down];
    kbs->right = to_depressed_state[kbs->right];
    kbs->left = to_depressed_state[kbs->left];
    kbs->a = to_depressed_state[kbs->a];
    kbs->s = to_depressed_state[kbs->s];
    kbs->d = to_depressed_state[kbs->d];
    kbs->w = to_depressed_state[kbs->w];
    kbs->q = to_depressed_state[kbs->q];
    kbs->e = to_depressed_state[kbs->e];
    kbs->r = to_depressed_state[kbs->r];

    kbs->one = to_depressed_state[kbs->one];
    kbs->two = to_depressed_state[kbs->two];
    kbs->three = to_depressed_state[kbs->three];
    kbs->four = to_depressed_state[kbs->four];
    kbs->five = to_depressed_state[kbs->five];
    kbs->six = to_depressed_state[kbs->six];
    kbs->seven = to_depressed_state[kbs->seven];
    kbs->eight = to_depressed_state[kbs->eight];
}

bool IsPressed(enum KeyboardKeyState state) {
    static bool lookup[KeyboardKeyState_COUNT] = {0, 0, 1, 1};
    return lookup[state];
}

bool IsNewlyPressed(enum KeyboardKeyState state) {
    return state == KeyboardKeyState_Pressed;
}
