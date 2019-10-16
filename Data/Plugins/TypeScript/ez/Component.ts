import __GameObject = require("./GameObject")
export import GameObject = __GameObject.GameObject;

declare function __CPP_Component_GetOwner(component: Component): GameObject;
declare function __CPP_Component_SetActive(component: Component, active: boolean): GameObject;

export abstract class Component
{
    GetOwner(): GameObject 
    {
        return __CPP_Component_GetOwner(this);
    }

    SetActive(active: boolean): void
    {
        __CPP_Component_SetActive(this, active);
    }
}

export abstract class TypescriptComponent extends Component
{
    abstract Update(): void;
}
