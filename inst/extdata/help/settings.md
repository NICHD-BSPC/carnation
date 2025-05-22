### Settings
-------------

The Settings tab provides administrative functionality for configuring Carnation's data sources and user access. This tab is only visible to users with administrative privileges in multi-user deployments.

#### Data Areas

Data areas are directories containing RNA-seq analyses that Carnation can access. Each data area is associated with a specific user group.

- **View Data Areas**: The table shows currently configured data areas and their associated user groups.
- **Add Data Area**: 
  - Click the "Add data area" button to configure a new data source.
  - Enter a user group name (existing or new) that should have access to this data.
  - Provide the full path to the directory containing your analyses.
  - The directory should follow Carnations's expected structure (see below).
  - Click "OK" to save the new data area.
- **Remove Data Area**: 
  - Select one or more data areas from the table.
  - Click the "Remove selected" button to remove them.
  - Confirm the deletion when prompted.

#### User Management

This section allows administrators to manage user access to different data areas. This is only shown in multi-user deployments.

- **View Users**: The table shows currently configured users and their associated groups.
- **Add User**:
  - Click the "Add user" button to add a new user.
  - Enter the username.
  - Select one or more user groups this user should belong to.
  - Click "OK" to save the new user.
- **Remove User**:
  - Select one or more users from the table.
  - Click the "Remove selected" button to remove them.
  - Confirm the deletion when prompted.

#### Data Directory Structure

Carnation expects a specific directory structure for data areas:

```
/carnation/data/
  ├─ project1
  │  ├─ main.rds
  │  ├─ subset.rds
  │  └─ main-nooutlier.rds
  │
  └─ project2
     └─ default.rds
```

#### Project Descriptions

You can add descriptions to projects by creating a `project-description.yaml` file in the project directory. For example, for the
above example, a `project-description.yaml` file in the `project1` directory could look like this:

```yaml
main: "Standard analysis of all samples"
subset: "Standard analysis using only samples passing QC"
main-nooutlier: "Standard analysis excluding outlier samples"
```

#### Access Control

- **User Groups**: Users can belong to multiple groups, giving them access to all data areas associated with those groups.
- **Admin Group**: Users in the admin group (configured in `config.yaml`) have access to the Settings tab and can manage data areas and users.
- **Changes**: After making changes to data areas or users, click "Save changes" to apply them. This will update the access configuration and may require a reload of the application.

#### Configuration Files

Carnation stores access configuration in an YAML file. The initial setup will prompt for:

1. A username (defaults to the current user)
2. A user group (defaults to "admin")
3. The path to the first data area

This configuration can be modified through the Settings tab after initial setup.

**Note**: Changes to data areas and user access require administrative privileges and will affect all users of the application.
