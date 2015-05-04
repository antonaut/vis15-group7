/****************************************************************************
** Meta object code from reading C++ file 'ObjectListElementEditor.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../ObjectListElementEditor.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ObjectListElementEditor.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ObjectListElementEditor[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      25,   24,   24,   24, 0x0a,
      40,   24,   24,   24, 0x0a,
      58,   24,   24,   24, 0x0a,
      70,   24,   24,   24, 0x0a,
      84,   24,   24,   24, 0x0a,
      94,   24,   24,   24, 0x0a,
     110,   24,   24,   24, 0x0a,
     122,   24,   24,   24, 0x0a,
     131,   24,   24,   24, 0x0a,
     156,  144,   24,   24, 0x0a,
     203,  199,   24,   24, 0x0a,
     232,   24,   24,   24, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_ObjectListElementEditor[] = {
    "ObjectListElementEditor\0\0sltMoveToTop()\0"
    "sltMoveToBottom()\0sltMoveUp()\0"
    "sltMoveDown()\0sltEdit()\0sltEditMember()\0"
    "sltDelete()\0sltAdd()\0sltAddAddr()\0"
    "item,column\0sltItemDoubleClicked(QTreeWidgetItem*,int)\0"
    "obj\0sltEditorDestroyed(QObject*)\0"
    "sltSubObjectModified()\0"
};

void ObjectListElementEditor::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ObjectListElementEditor *_t = static_cast<ObjectListElementEditor *>(_o);
        switch (_id) {
        case 0: _t->sltMoveToTop(); break;
        case 1: _t->sltMoveToBottom(); break;
        case 2: _t->sltMoveUp(); break;
        case 3: _t->sltMoveDown(); break;
        case 4: _t->sltEdit(); break;
        case 5: _t->sltEditMember(); break;
        case 6: _t->sltDelete(); break;
        case 7: _t->sltAdd(); break;
        case 8: _t->sltAddAddr(); break;
        case 9: _t->sltItemDoubleClicked((*reinterpret_cast< QTreeWidgetItem*(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 10: _t->sltEditorDestroyed((*reinterpret_cast< QObject*(*)>(_a[1]))); break;
        case 11: _t->sltSubObjectModified(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData ObjectListElementEditor::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ObjectListElementEditor::staticMetaObject = {
    { &ElementEditor::staticMetaObject, qt_meta_stringdata_ObjectListElementEditor,
      qt_meta_data_ObjectListElementEditor, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ObjectListElementEditor::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ObjectListElementEditor::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ObjectListElementEditor::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ObjectListElementEditor))
        return static_cast<void*>(const_cast< ObjectListElementEditor*>(this));
    return ElementEditor::qt_metacast(_clname);
}

int ObjectListElementEditor::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = ElementEditor::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
