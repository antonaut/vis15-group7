/****************************************************************************
** Meta object code from reading C++ file 'ObjectElementEditor.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../ObjectElementEditor.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ObjectElementEditor.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ObjectElementEditor[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      21,   20,   20,   20, 0x0a,
      37,   20,   20,   20, 0x0a,
      66,   59,   20,   20, 0x0a,
      94,   90,   20,   20, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_ObjectElementEditor[] = {
    "ObjectElementEditor\0\0sltEditButton()\0"
    "sltViewObjectButton()\0action\0"
    "sltHandleMenu(QAction*)\0obj\0"
    "sltDestroyEditorWidget(QObject*)\0"
};

void ObjectElementEditor::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ObjectElementEditor *_t = static_cast<ObjectElementEditor *>(_o);
        switch (_id) {
        case 0: _t->sltEditButton(); break;
        case 1: _t->sltViewObjectButton(); break;
        case 2: _t->sltHandleMenu((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 3: _t->sltDestroyEditorWidget((*reinterpret_cast< QObject*(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData ObjectElementEditor::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ObjectElementEditor::staticMetaObject = {
    { &ElementEditor::staticMetaObject, qt_meta_stringdata_ObjectElementEditor,
      qt_meta_data_ObjectElementEditor, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ObjectElementEditor::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ObjectElementEditor::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ObjectElementEditor::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ObjectElementEditor))
        return static_cast<void*>(const_cast< ObjectElementEditor*>(this));
    return ElementEditor::qt_metacast(_clname);
}

int ObjectElementEditor::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = ElementEditor::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
