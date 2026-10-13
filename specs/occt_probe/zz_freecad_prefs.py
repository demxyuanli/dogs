import sys

import FreeCAD

def dump(path):
    p = FreeCAD.ParamGet(path)
    print("=== %s ===" % path)
    try:
        for k in p.GetBools():
            print("  bool  %s = %s" % (k, p.GetBool(k)))
        for k in p.GetInts():
            print("  int   %s = %s" % (k, p.GetInt(k)))
        for k in p.GetFloats():
            print("  float %s = %s" % (k, p.GetFloat(k)))
        for k in p.GetStrings():
            print("  str   %s = %s" % (k, p.GetString(k)))
    except Exception as e:
        print("  err", e)

dump("User parameter:BaseApp/Preferences/Mod/Part")
dump("User parameter:BaseApp/Preferences/Mod/Part/ImportIges")
dump("User parameter:BaseApp/Preferences/Mod/Part/ImportStep")
sys.stdout.flush()
