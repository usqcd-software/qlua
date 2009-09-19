#include <qlua.h>                                                    /* DEPS */
#include <qxml.h>                                                    /* DEPS */
#include <string.h>

static const char xml_lib[] = "xml";

typedef const char *(parse_func)(lua_State *, const char*, const char*, char *);

/* predicates */
static int
is_doctypedecl(const char *s, const char *e)
{
    return ((e > s + 9) && (strncmp(s, "<!DOCTYPE", 9) == 0));
}

static int
is_Comment(const char *s, const char *e)
{
    return ((e > s + 4) && (strncmp(s, "<!--", 4) == 0));
}

static int
is_PI(const char *s, const char *e)
{
    return ((e > s + 2) && (strncmp(s, "<?", 2) == 0));
} 

static int
is_S(const char *s, const char *e)
{
    if (e > s)
        return ((*s == ' ') || (*s == '\t') || (*s == '\n') || (*s == '\r'));
    else
        return 0;
}

static int
is_CDSect(const char *s, const char *e)
{
    return ((e > s + 9) && (strncmp(s, "<![CDATA[", 9) == 0));
}

static int
is_ETag(const char *s, const char *e)
{
    return ((e > s + 2) && (strncmp(s, "</", 2) == 0));
}

static int
is_XMLDecl(const char *start, const char *end)
{
    return ((end > start + 5) && (strncmp(start, "<?xml", 5) == 0));
}

static int
is_Misc(const char *start, const char *end)
{
    return (is_Comment(start, end) ||
            is_PI(start, end) ||
            is_S(start, end));
}

/* XML rule parsers */
static const char *
parse_S(lua_State *L,
        const char *s,
        const char *e,
        char *buf)
{
    size_t sp = strspn(s, " \t\n\r");
    if (sp == 0)
        luaL_error(L, "xml.parse() expecting S");

    return s + sp;
}

static const char *
parse_S_opt(lua_State *L,
            const char *s,
            const char *e,
            char *buf)
{
    if (is_S(s, e))
        return parse_S(L, s, e, buf);
    return s;
}

static const char *
parse_Comment(lua_State *L,
              const char *s,
              const char *e,
              char *buf)
{
    static const char err[] = "xml.parse() expecting Comment";
    const char *pos;

    if ((e < s + 7) || (strncmp(s, "<!--", 4) != 0))
        luaL_error(L, err);

    for (pos = s + 4; pos < e; pos++) {
        if (*pos != '-')
            continue;
        if (pos + 3 > e)
            luaL_error(L, err);
        if (pos[1] != '-')
            continue;
        if (pos[2] != '>')
            luaL_error(L, err);
        return pos + 3;
    }
    luaL_error(L, err);
    return NULL; /* can not happen */
}

static const char *
parse_Name(lua_State *L,
           const char *s,
           const char *e,
           char *buf)
{
    const char *pos = s;

    if ((s >= e) ||
        !((*s == ':') ||
          (*s == '_') ||
          ((*s >= 'A') && (*s <= 'Z')) ||
          ((*s >= 'a') && (*s <= 'z'))))
        luaL_error(L, "xml.parse() expecting Name");
    *buf++ = *s;
    for (pos = s + 1; pos < e; pos++, buf++) {
        if (!((*pos == ':') ||
              (*pos == '_') ||
              ((*pos >= 'A') && (*pos <= 'Z')) ||
              ((*pos >= 'a') && (*pos <= 'z')) ||
              ((*pos >= '0') && (*pos <= '9'))))
            break;
        *buf = *pos;
    }
    *buf = 0;

    return pos;
}

static const char *
parse_PITarget(lua_State *L,
               const char *s,
               const char *e,
               char *buf)
{
    const char *pos = parse_Name(L, s, e, buf);
    
    if (strcasecmp(buf, "xml") == 0)
        luaL_error(L, "xml.parse() expecting PITarget");

    return pos;
}

static const char *
parse_PI(lua_State *L,
         const char *s,
         const char *e,
         char *buf)
{
    static const char err[] = "xml.parse() expecting PI";
    const char *pos;

    if ((e < s + 4) || (strncmp(s, "<?", 2) != 0))
        luaL_error(L, err);

    pos = parse_PITarget(L, s + 2, e, buf);
    if (is_S(pos, e)) {
        for (; pos < e; pos++) {
            if (*pos != '?')
                continue;
            if (pos + 2 > e)
                luaL_error(L, err);
            if (pos[1] != '>')
                continue;
            return pos + 2;
        }
        luaL_error(L, err);
        return NULL; /* never happens */
    } else {
        if ((pos + 2 > e) || (strncmp(pos, "?>", 2) != 0))
            luaL_error(L, err);
        return pos + 2;
    }
}

static const char *
parse_VersionNum(lua_State *L,
                 const char *s,
                 const char *e,
                 char *buf)
{
    const char *pos;

    if ((s + 2>= e) || (strncmp(s, "1.", 2) != 0))
        luaL_error(L, "xml.parse() expecting VersionNum");
    for (pos = s + 2; pos < e; pos++) {
        if ((*pos < '0') || (*pos > '9'))
            break;
    }

    return pos;
}

static const char *
parse_Eq(lua_State *L,
         const char *s,
         const char *e,
         char *buf)
{
    const char *pos = s;

    pos = parse_S_opt(L, pos, e, buf);
    if ((pos >= e) || (*pos != '='))
        luaL_error(L, "xml.parse() expecting Eq");
    pos = parse_S_opt(L, pos + 1, e, buf);

    return pos;
}

static const char *
parse_VersionInfo(lua_State *L,
                  const char *s,
                  const char *e,
                  char *buf)
{
    const char err[] = "xml.parse() expecting VersionInfo";
    const char *pos;

    pos = parse_S(L, s, e, buf);
    if ((pos + 7 > e) || (strncmp(pos, "version", 7) != 0))
        luaL_error(L, err);
    pos = parse_Eq(L, pos + 7, e, buf);
    if (pos >= e)
        luaL_error(L, err);
    switch (*pos) {
    case '\'':
        pos = parse_VersionNum(L, pos + 1, e, buf);
        if ((pos + 1 >= e) || (*pos != '\''))
            luaL_error(L, err);
        return pos + 1;
    case '\"':
        pos = parse_VersionNum(L, pos + 1, e, buf);
        if ((pos + 1 >= e) || (*pos != '\"'))
            luaL_error(L, err);
        return pos + 1;
    default:
        luaL_error(L, err);
        return NULL; /* never happens */
    }
    
}

static const char *
parse_XMLDecl(lua_State *L,
              const char *s,
              const char *e,
              char *buf)
{
    static const char err[] = "xml.parse() expecting XMLDecl";
    const char *pos;
    
    if ((e < s + 8) || strncmp(s, "<?xml", 5) != 0)
        luaL_error(L, err);
    pos = parse_VersionInfo(L, s + 5, e, buf);
    pos = parse_S_opt(L, pos, e, buf);
    if ((e < pos + 2) || strncmp(pos, "?>", 2) != 0)
        luaL_error(L, err);

    return pos + 2;
}

static const char *
parse_Misc(lua_State *L,
           const char *s,
           const char *e,
           char *buf)
{
    if (is_Comment(s, e))
        return parse_Comment(L, s, e, buf);
    if (is_PI(s, e))
        return parse_PI(L, s, e, buf);
    if (is_S(s, e))
        return parse_S(L, s, e, buf);

    luaL_error(L, "xml.parse(): illformed Misc");
    return NULL; /* can not happen */
}

static const char *
parse_Misc_opt(lua_State *L,
               const char *s,
               const char *e,
               char *buf)
{
    while (is_Misc(s, e))
        s = parse_Misc(L, s, e, buf);

    return s;
}

static void
check_valid_char(lua_State *L, unsigned int v)
{
    if (v >= 0x20 && v <= 0x7f)
        return;
    if (v == 0x9 || v == 0xa || v == 0xd)
        return;
    luaL_error(L, "xml.parse() character out of repertoire");
}

static const char *
parse_ETag(lua_State *L,
           const char *s,
           const char *e,
           char *buf)
{
    const char err[] = "xml.parse() expecting ETag";
    const char *pos;
    const char *tag = luaL_checkstring(L, -1);

    if (!is_ETag(s, e))
        luaL_error(L, err);
    pos = parse_Name(L, s + 2, e, buf);
    if (strcmp(tag, buf) != 0)
        luaL_error(L, err);
    lua_pop(L, 1);
    pos = parse_S_opt(L, pos, e, buf);
    if ((pos >= e) || (*pos != '>'))
        luaL_error(L, err);

    return pos + 1;
}

static const char *
parse_CharRef(lua_State *L,
              const char *s,
              const char *e,
              char *buf)
{
    const char err[] = "xml.parse() expecting CharRef";
    int f = 1;
    unsigned int v = 0;

    if (s >= e)
        luaL_error(L, err);
    switch (*s) {
    case 'x':
        for (s++ ; s < e; s++) {
            if (*s == ';') {
                if (f)
                    luaL_error(L, err);
                check_valid_char(L, v);
                *buf = v;
                return s + 1;
            } else if (*s >= '0' && *s <= '9') {
                v = v * 16 + *s - '0';
            } else if (*s >= 'a' && *s <= 'f') {
                v = v * 16 + *s + 10 - 'a';
            } else if (*s >= 'A' && *s <= 'F') {
                v = v * 16 + *s + 10 - 'A';
            } else {
                luaL_error(L, err);
            }
            if (v > 0x7f)
                luaL_error(L, err);
            f = 0;
        }
        break;
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9':
        v = *s - '0';
        for (s++; s < e; s++) {
            if (*s == ';') {
                check_valid_char(L, v);
                *buf = v;
                return s + 1;
            } else if (*s >= '0' && *s <= '9') {
                v = v * 10 + *s - '0';
            } else {
                luaL_error(L, err);
            }
            if (v > 0x7f)
                luaL_error(L, err);
        }
        break;
    }
    luaL_error(L, err);
    return NULL; /* never happens */
}

static const char *
parse_Reference(lua_State *L,
                const char *s,
                const char *e,
                char *buf)
{
    const char err[] = "xml.parse() expecting Reference";
    const char *pos = s;
    char v = 0;

    if (s >= e)
        luaL_error(L, err);
    if (*s == '#')
        return parse_CharRef(L, s + 1, e, buf);
    if ((s + 5 < e) && (strncmp(s, "quot;", 5) == 0))
        v = '\"', pos = s + 5;
    else if ((s + 5 < e) && (strncmp(s, "apos;", 5) == 0))
        v = '\'', pos = s + 5;
    else if ((s + 4 < e) && (strncmp(s, "amp;", 4) == 0))
        v = '&', pos = s + 4;
    else if ((s + 3 < e) && (strncmp(s, "lt;", 3) == 0))
        v = '<', pos = s + 3;
    else if ((s + 3 < e) && (strncmp(s, "gt;", 3) == 0))
        v = '>', pos = s + 3;
    else
        luaL_error(L, err);
    *buf = v;

    return pos;
}

static const char *
parse_AttValue(lua_State *L,
              const char *s,
              const char *e,
              char *buf)
{
    const char err[] = "xml.parse() expecting AttValue";
    const char *pos = s;

    if (s >= e)
        luaL_error(L, err);
    switch (*s) {
    case '\'':
    case '\"':
        for (pos = s + 1; pos < e; pos++) {
            if (*pos == *s) {
                *buf = 0;
                return pos + 1;
            }
            switch (*pos) {
            case '<':
                luaL_error(L, err);
                break;
            case '&':
                pos = parse_Reference(L, pos + 1, e, buf) - 1;
                buf++;
                break;
            default:
                check_valid_char(L, *pos);
                *buf++ = *pos;
                break;
            }
        }
        /* through */
    default:
        luaL_error(L, err);
        break;
    }
    return NULL; /* can not happen */
}

static const char *
parse_Attribute(lua_State *L,
                const char *s,
                const char *e,
                char *buf)
{
    const char *pos;

    pos = parse_Name(L, s, e, buf);
    lua_pushstring(L, buf);
    pos = parse_Eq(L, pos, e, buf);
    pos = parse_AttValue(L, pos, e, buf);
    lua_pushstring(L, buf);
    lua_setfield(L, -4, luaL_checkstring(L, -2));
    lua_pop(L, 1);

    return pos;
}

static void
add_chardata_segment(lua_State *L, const char *buf)
{
    int len = lua_objlen(L, -2);

    if (*buf  == 0)
        return;
    if (len == 0) {
        lua_pushstring(L, buf);
        lua_rawseti(L, -3, 1);
    } else {
        lua_rawgeti(L, -2, len);
        if (lua_isstring(L, -1)) {
            lua_pushstring(L, buf);
            lua_concat(L, 2);
            lua_rawseti(L, -3, len);
        } else {
            lua_pop(L, 1);
            lua_pushstring(L, buf);
            lua_rawseti(L, -3, len + 1);
        }
    }
}

static const char *
parse_CharData_opt(lua_State *L,
                   const char *s,
                   const char *e,
                   char *buf)
{
    char *bs = buf;

    for (; s < e; s++, buf++) {
        if (*s == '<')
            break;
        if (*s == '&')
            break;
        if ((*s == ']') && (s + 3 < e) && (strncmp(s, "]]>", 3) == 0))
            luaL_error(L, "xml.parse() expecting CharData");
        check_valid_char(L, *s);
        *buf = *s;
    }
    *buf = 0;
    add_chardata_segment(L, bs);
    return s;
}

static const char *
parse_CDSect(lua_State *L,
             const char *s,
             const char *e,
             char *buf)
{
    char *bs = buf;
    
    for (s += 9; s < e; s++, buf++) {
        if ((*s == ']') && (s + 3 <= e) && (strncmp(s, "]]>", 3) == 0)) {
            *buf = 0;
            add_chardata_segment(L, bs);
            return s + 3;
        }
        check_valid_char(L, *s);
        *buf = *s;
    }
    luaL_error(L, "xml.parse() expecting CDSect");

    return NULL;
}

static parse_func parse_element;

static const char *
parse_content(lua_State *L, /* [-1] = tag, [-2] = table */
              const char *s,
              const char *e,
              char *buf)
{
    const char *pos = parse_CharData_opt(L, s, e, buf);

    while (pos < e) {
        switch (*pos) {
        case '<':
            if (is_Comment(pos, e))
                pos = parse_Comment(L, pos, e, buf);
            else if (is_PI(pos, e))
                pos = parse_PI(L, pos, e, buf);
            else if (is_CDSect(pos, e))
                pos = parse_CDSect(L, pos, e, buf);
            else if (is_ETag(pos, e))
                return pos;
            else {
                pos = parse_element(L, pos, e, buf);
                lua_rawseti(L, -3, lua_objlen(L, -3) + 1);
            }
            pos = parse_CharData_opt(L, pos, e, buf);
            break;
        case '&':
            pos = parse_Reference(L, pos + 1, e, buf);
            buf[1] = 0;
            add_chardata_segment(L, buf);
            pos = parse_CharData_opt(L, pos, e, buf);
            break;
        default:
            goto error;
        }
        
    }
error:
    luaL_error(L, "xml.parse() expecting content");
    return NULL;
}

static const char *
parse_element(lua_State *L,
              const char *s,
              const char *e,
              char *buf)
{
    const char err[] = "xml.parse() expecting element";
    const char *pos;

    if ((s >= e) || (*s != '<'))
        luaL_error(L, err);

    pos = parse_Name(L, s + 1, e, buf);
    lua_newtable(L);
    lua_pushstring(L, buf);
    lua_setfield(L, -2, "$tag");
    lua_pushstring(L, buf);
    for (;;) {
        pos = parse_S_opt(L, pos, e, buf);
        if (pos >= e)
            luaL_error(L, err);
        switch (*pos) {
        case '>':
            pos = parse_content(L, pos + 1, e, buf);
            return parse_ETag(L, pos, e, buf);
        case '/':
            if ((pos + 1 >= e) || (pos[1] != '>'))
                luaL_error(L, err);
            lua_pop(L, 1);
            return pos + 2;
        default:
            pos = parse_Attribute(L, pos, e, buf);
            break;
        }
    }
}

static const char *
parse_prolog(lua_State *L,
             const char *s,
             const char *e,
             char *buf)
{
    const char *pos = s;

    if (is_XMLDecl(pos, e))
        pos = parse_XMLDecl(L, pos, e, buf);
    while (is_Misc(pos, e))
        pos = parse_Misc(L, pos, e, buf);
    if (is_doctypedecl(pos, e))
        luaL_error(L, "xml.parse() doctypedecl is not supported");

    return pos;
}

static void
parse_document(lua_State *L,
               const char *s, 
               const char *e,
               char *buf)
{
    const char *pos;

    pos = parse_prolog(L, s, e, buf);
    pos = parse_element(L, pos, e, buf);
    pos = parse_Misc_opt(L, pos, e, buf);
    if (pos != e)
        luaL_error(L, "illformed XML document");
}

static int
q_xml_parse(lua_State *L)
{
    const char *in = luaL_checkstring(L, 1);
    const size_t len = lua_objlen(L, 1);
    char *buf = qlua_malloc(L, len);

    parse_document(L, in, in + len, buf);
    qlua_free(L, buf);
    
    return 1;
}

static int
q_xml_unparse(lua_State *L)
{
    return luaL_error(L, "XXX xml.unparse() not implemented");
}

static struct luaL_Reg fXML[] = {
    { "parse",    q_xml_parse },
    { "unparse",  q_xml_unparse },
    { NULL,       NULL }
};

int
init_xml(lua_State *L)
{
    luaL_register(L, xml_lib, fXML);
    return 0;
}

int
fini_xml(lua_State *L)
{
    return  0;
}
